"""Scientific failure modes and complete workflow contracts."""

import csv
import json
import shutil
import time
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from proteinwork.diagnostics import diagnose, summarize_observable
from proteinwork.direct import SamplingResult, sample_direct
from proteinwork.ingest import ASSAY_FIELDS, ingest_measurements
from proteinwork.model import PottsModel, restrict_model
from proteinwork.planner import freeze_plan
from proteinwork.validation import apply_residue_requests, digest, load_request, read_json, save_json
from proteinwork.workflow import _runtime, _sample_with_checkpoints, run_workflow, verify_run


ROOT = Path(__file__).resolve().parents[1]


def request():
    return {"schema_version": 1, "benchmark_config": "benchmark.json", "evaluation_config": "evaluation.json",
            "case_id": "cm_dev_exact_8_12", "sampling": {"method": "direct", "seeds": [1, 2, 3, 4],
            "draws": 2048, "warmup": 32}, "selection": {"batch_size": 2, "max_pool": 4, "optimizer": "milp"},
            "allow_unqualified_export": True}


def fixture():
    model = PottsModel("AC", np.array([[0.0, 0.2], [0.1, -0.1]]), np.zeros((2, 2, 2, 2)))
    config = {"archive": {"member_sha256": "fixture_model"}, "model": {"beta": 1.0},
              "coordinate_reference": {"uniprot_accession": "fixture", "reference_span_in_full_accession": [2, 3]}}
    mapping = [{"model_index_0based": i, "alignment_column_1based": i + 1, "anchor_symbol": "A",
                "uniprot_position_1based": i + 2} for i in range(2)]
    return SimpleNamespace(model=model, anchor="AA", witness="AA", mapping=mapping, config=config,
        cases=[{"case_id": "cm_dev_exact_8_12", "kind": "exact_conditional", "split": "development"}],
        provenance={"fixture": True})


class DiagnosticTests(unittest.TestCase):
    def test_checkpoints_continue_states_and_only_warm_up_once(self):
        inputs = fixture()
        reduced = restrict_model(inputs.model, [(0, 1), (0, 1)])
        evaluation = read_json(ROOT / "configs/evaluation.json")
        evaluation["timing"]["checkpoints_post_warmup_observations_per_chain"] = [32, 64, 96]
        sampling = {"method": "direct", "seeds": [10, 20, 30, 40], "draws": 96,
                    "warmup": 8, "checkpoint_policy": "first_qualified"}
        diagnoses = [{"status": "accuracy_not_reached", "checks": {"fixture": False}},
                     {"status": "accuracy_passed", "checks": {"fixture": True}}]
        with patch("proteinwork.workflow.diagnose", side_effect=diagnoses), patch(
            "proteinwork.workflow.sample_direct", wraps=sample_direct) as sampler:
            result, diagnosis, _, _, records = _sample_with_checkpoints(reduced, sampling=sampling, beta=1,
                anchor="AA", selected_pairs=[], evaluation=evaluation, reference=None, started=time.perf_counter(),
                budget_seconds=100)
        self.assertEqual(result.states.shape, (4, 64, 2))
        self.assertEqual([r["draws_per_chain"] for r in records], [32, 64])
        self.assertEqual(sampler.call_args_list[0].kwargs["warmup"], 8)
        self.assertEqual(sampler.call_args_list[1].kwargs["warmup"], 0)
        np.testing.assert_array_equal(sampler.call_args_list[1].kwargs["initial_states"], result.states[:, 31])
        self.assertFalse(set(records[0]["segment_seeds"]) & set(records[1]["segment_seeds"]))

    def test_independent_and_stuck_traces_are_distinguished(self):
        rng = np.random.default_rng(77)
        normal = summarize_observable(rng.normal(size=(4, 4000)))
        self.assertLess(normal["rhat"], 1.01)
        self.assertGreater(normal["ess_bulk"], 400)
        self.assertGreater(normal["ess_tail"], 400)
        stuck = summarize_observable(np.ones((4, 1000)))
        self.assertEqual(stuck["status"], "insufficient_or_constant_trace")
        self.assertFalse(stuck["constant_proven"])
        self.assertIsNone(stuck["mcse"])

    def test_known_constant_target_and_wrong_value(self):
        self.assertEqual(summarize_observable(np.ones((4, 20)), proven_value=1)["status"], "constant_verified")
        self.assertEqual(summarize_observable(np.ones((4, 20)), proven_value=0)["status"], "constant_value_mismatch")

    def test_zero_variable_target_is_qualified_without_invented_ess(self):
        inputs = fixture()
        model = restrict_model(inputs.model, [(0,), (0,)])
        result = sample_direct(model, beta=1, seeds=[1, 2, 3, 4], draws=16, warmup=0)
        evaluation = read_json(ROOT / "configs/evaluation.json")
        diagnosis = diagnose(model, result, beta=1, anchor="AA", selected_pairs=[], evaluation=evaluation)
        self.assertEqual(diagnosis["status"], "accuracy_passed")
        self.assertTrue(diagnosis["observables"]["energy"]["constant_proven"])
        self.assertIsNone(diagnosis["observables"]["energy"]["ess_bulk"])
        self.assertEqual(diagnosis["comparison"]["total_variation"], 0)

    def test_stuck_sample_does_not_pass_exact_accuracy(self):
        inputs = fixture()
        model = restrict_model(inputs.model, [(0, 1), (0, 1)])
        states = np.zeros((4, 1024, 2), dtype=int)
        result = SamplingResult(states, model.energy(states), states[:, 0], np.zeros((4, 1024), bool), {}, {})
        diagnosis = diagnose(model, result, beta=1, anchor="AA", selected_pairs=[[0, 1]],
                             evaluation=read_json(ROOT / "configs/evaluation.json"))
        self.assertEqual(diagnosis["status"], "accuracy_not_reached")
        self.assertFalse(diagnosis["checks"]["energy_diagnostics"])


class InputRulesTests(unittest.TestCase):
    def test_windows_processor_identity_survives_empty_platform_probe(self):
        with patch.dict("os.environ", {"PROCESSOR_IDENTIFIER": "fixture_cpu"}), patch(
            "proteinwork.workflow.platform.processor", return_value=""):
            self.assertEqual(_runtime()["processor"], "fixture_cpu")

    def test_duplicate_unknown_and_ambiguous_json_fail(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "request.json"
            for text in ('{"a":1,"a":2}', '{"a": NaN}', '[]'):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_json(path)
            for key, value in (("execution", "physical"), ("unknown_rule", True)):
                save_json(path, {**request(), key: value})
                with self.assertRaises(ValueError):
                    load_request(path)
            invalid = request()
            invalid["sampling"]["seeds"] = [True, 2, 3, 4]
            save_json(path, invalid)
            with self.assertRaises(ValueError):
                load_request(path)

    def test_numbering_conflicts_and_forbidden_requests_fail(self):
        inputs = fixture()
        domains = ((0,), (0, 1))
        good = {"uniprot_position_1based": 3, "accession": "fixture", "allowed_symbols": "C",
                "expected_anchor_symbol": "A", "reason": "Test restriction"}
        result, _ = apply_residue_requests(inputs, domains, [good])
        self.assertEqual(result, ((0,), (1,)))
        for bad in ([good, good], [{**good, "accession": "other"}],
                    [{**good, "uniprot_position_1based": 2}], [{**good, "expected_anchor_symbol": "C"}]):
            with self.assertRaises(ValueError):
                apply_residue_requests(inputs, domains, bad)


class WorkflowTests(unittest.TestCase):
    def _run(self, directory, allow=True):
        inputs = fixture()
        cfg = request()
        cfg["allow_unqualified_export"] = allow
        path = directory / "request.json"
        save_json(path, cfg)
        evaluation = read_json(ROOT / "configs/evaluation.json")
        evaluation["selected_pair_observables"][cfg["case_id"]] = [[0, 1]]
        save_json(directory / "evaluation.json", evaluation)
        save_json(directory / "benchmark.json", inputs.config)
        with patch("proteinwork.workflow.load_benchmark", return_value=inputs), patch(
            "proteinwork.workflow.benchmark_domains", return_value=((0, 1), (0, 1))):
            manifest = run_workflow(path, directory / "run")
        return manifest

    def test_complete_outputs_repeats_identity_and_assay_round_trip(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest = self._run(directory)
            run = directory / "run"
            self.assertEqual(manifest["raw_observations"], 8192)
            self.assertEqual(manifest["exported_candidates"], 2)
            self.assertEqual(verify_run(run)["status"], "verified")
            with (run / "assay_worksheet.csv").open(newline="") as handle:
                reader = csv.DictReader(handle)
                columns, rows = reader.fieldnames, list(reader)
            for index, row in enumerate(rows):
                row.update(value=str(index + 1), assay_id="test_fixture", units="arbitrary", context="synthetic_unit_test", round="0")
            measurements = directory / "measurements.csv"
            with measurements.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns)
                writer.writeheader()
                writer.writerows(rows)
            receipt = ingest_measurements(run, measurements, directory / "assays.csv")
            self.assertEqual(receipt["new_measurements"], 2)
            with (directory / "assays.csv").open(newline="") as handle:
                self.assertEqual(csv.DictReader(handle).fieldnames, ASSAY_FIELDS)
            with (run / "candidates.fasta").open("a") as handle:
                handle.write("corruption")
            with self.assertRaises(ValueError):
                verify_run(run)

    def test_unqualified_export_is_blocked_and_existing_run_is_preserved(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp)
            with patch("proteinwork.workflow.diagnose", return_value={"status": "accuracy_not_reached", "checks": {}}):
                manifest = self._run(directory, allow=False)
            self.assertEqual(manifest["exported_candidates"], 0)
            self.assertEqual((directory / "run/candidates.fasta").read_text(), "")
            with self.assertRaises(ValueError):
                run_workflow(directory / "request.json", directory / "run")

    def test_evaluation_evidence_cannot_train_planner(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest = self._run(directory)
            manifest["case"]["split"] = "evaluation"
            save_json(directory / "run/manifest.json", manifest)
            with self.assertRaises(ValueError):
                freeze_plan([directory / "run"], directory / "plan.json")

    def test_planner_keeps_methods_separate_and_rejects_duplicate_evidence(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest = self._run(directory)
            manifest["elapsed_seconds"] = 10.0
            save_json(directory / "run/manifest.json", manifest)
            second = directory / "second"
            shutil.copytree(directory / "run", second)
            alternate = dict(manifest)
            alternate["method"], alternate["elapsed_seconds"] = "encoded", 20.0
            save_json(second / "manifest.json", alternate)
            plan = freeze_plan([directory / "run", second], directory / "plan.json")
            self.assertEqual(plan["choices"]["exact_conditional"]["method"], "direct")
            self.assertEqual(plan["choices"]["exact_conditional"]["median_seconds"], 10.0)
            with self.assertRaises(ValueError):
                freeze_plan([directory / "run", directory / "run"], directory / "duplicate.json")
            alternate["comparison_context_sha256"] = "different_hardware_or_workload"
            save_json(second / "manifest.json", alternate)
            with self.assertRaises(ValueError):
                freeze_plan([directory / "run", second], directory / "mismatched.json")


if __name__ == "__main__":
    unittest.main()

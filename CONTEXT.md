# Protein candidate discovery

**Research objective**: The protein property a researcher wants to investigate, together with its measurement, units and experimental context. A family model score is a computational preference and cannot substitute for this measurement.

**Coordinate anchor**: The unchanged reference sequence that defines residue numbering and sequence distance. An anchor may lie outside the trained model's support.

**Supported sequence**: A sequence whose selected residues all have finite original model fields and satisfy the declared sampling restrictions.

**Sampling target**: The conditional probability distribution defined by the immutable protein model, allowed residues and computational temperature.

**Raw trace**: Every recorded chain state, including repeated states after rejection. Removing repeats changes the empirical sampling distribution.

**Candidate pool**: Distinct sequences extracted from a trace, with links to their original chain observations. Its frequencies are not estimates of the sampling target.

**Candidate library**: A selected group of sequences satisfying the researcher's selection rules and experimental batch budget. Selection changes the distribution represented by the raw trace.

**Assay record**: A measurement attached to an exact sequence, experimental context, readout and units. Similar sequences do not inherit its measured value.

**Encoded proposal**: A proposed sequence change produced by a binary representation, corrected against the original sampling target before retention in the trace.

**Physical execution**: A run on an identified Extropic device through a verified runtime. THRML simulation and an exported Ising problem do not establish physical execution.

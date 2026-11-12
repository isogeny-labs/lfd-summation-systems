import ICFallDegrees.SparsePrimeCertificate

open ICFallDegrees SparsePolynomial SparseCertificate SparsePrimeCertificate

-- Executable diagnostics only; theorem acceptance is independently checked
-- by the kernel in SparsePrimeCertificate, never inferred from this output.
#eval [
  (check normalizedInputs 3 normalizedUnit).map degree,
  (check rawInputs 4 rawUnit).map degree,
  (check normalizedInputs 2 normalizedUnit).map degree,
  (check rawInputs 3 rawUnit).map degree]

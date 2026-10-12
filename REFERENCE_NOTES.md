# Reference sequence review notes

- `markers.fna` contains 20 FASTA records supplied for this project.
- The dash in the supplied `espK` record was removed as a presumed alignment gap to make the query an ungapped DNA sequence. Verify this against the source sequence before relying on results.
- The standalone `A` line inside `Z2098` remains in the sequence as normal FASTA sequence-line wrapping. Verify the expected sequence length and source accession.
- Add accession, source publication/database, retrieval date/version, and redistribution terms for each reference before public release.
- Intimin reference hits are screened as separate marker columns. The pipeline does not declare an intimin subtype based on these hits alone.

- LOSAT Web's BLASTX takes a query and a subject subrange (`-query_loc` / `-subject_loc`), as the other modes do, and
  gives the bytes of the command line: each record is cut to its interval, a query record whose interval starts more
  than one letter past its end is not searched, and the results show the coordinates and frames of the whole
  records.

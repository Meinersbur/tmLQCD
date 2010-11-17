#include <config.h>
#include <lemon.h>
#include <memory.h>

int lemonReaderSetState(LemonReader *rdest, LemonReader *rsrc)
{
  MPI_Offset   disp;
  MPI_Datatype etype;
  MPI_Datatype ftype;
  char         drep[32];

  /* Set rdest reader state from rsrc */
  /* We do not copy the file pointer member fp - this needs to be set at construction */
  if (rdest->curr_header == (LemonRecordHeader*)NULL)
    rdest->curr_header = (LemonRecordHeader*)malloc(sizeof(LemonRecordHeader));
  if (!rsrc->is_awaiting_header)
    memcpy(rdest->curr_header, rsrc->curr_header, sizeof(LemonRecordHeader));

  rdest->is_last            = rsrc->is_last;
  rdest->is_awaiting_header = rsrc->is_awaiting_header;
  rdest->is_busy            = 0;
  rdest->is_striped         = 0;

  rdest->off                = rsrc->off;
  rdest->pos                = rsrc->pos;

  /* Now make the system agree with the reader state */
  MPI_File_get_view(*rsrc->fp, &disp, &etype, &ftype, drep);
  MPI_File_set_view(*rdest->fp, disp, etype, ftype, drep, MPI_INFO_NULL);
  MPI_File_seek(*rdest->fp, rdest->pos, MPI_SEEK_SET);

  MPI_Comm_dup(rsrc->cartesian, &rdest->cartesian);
  MPI_Comm_rank(rdest->cartesian, &rdest->my_rank);

  return LEMON_SUCCESS;
}

#include <config.h>
#include <lemon.h>
#include <stdio.h>

int lemonWriterSeek(LemonWriter *writer, MPI_Offset offset, int whence)
{
  int err;

  if (writer == (LemonWriter*)NULL || writer->is_awaiting_header)
  {
    fprintf(stderr, "[LEMON] Node %d reports in lemonWriterSeek:\n"
                    "        NULL pointer or uninitialized writer provided.\n", writer->my_rank);
    return LEMON_ERR_PARAM;
  }

  if (writer->is_busy)
    lemonFinishWriting(writer);

  if (whence == MPI_SEEK_CUR)
    writer->pos += offset;
  else if (whence == MPI_SEEK_SET)
    writer->pos = offset;
  else if (whence == MPI_SEEK_END)
    writer->pos = writer->data_length - offset;
  else
    return LEMON_ERR_PARAM;

  err = MPI_File_seek(*writer->fp, writer->pos, MPI_SEEK_SET);
  MPI_Barrier(writer->cartesian);

  if (err != MPI_SUCCESS)
  {
    fprintf(stderr, "[LEMON] Node %d reports in lemonWriterSeek:\n"
                    "        MPI_File_seek returned error %d.\n", writer->my_rank, err);
    return LEMON_ERR_SEEK;
  }
  return LEMON_SUCCESS;
}

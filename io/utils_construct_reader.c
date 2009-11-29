#include "utils.ih"

void construct_reader(READER ** reader, char const *filename)
{
  LIME_FILE *fh = NULL;
  int status = 0;

#ifdef HAVE_LIBLEMON
  fh = (MPI_File*)malloc(sizeof(MPI_File));
  status = MPI_File_open(g_cart_grid, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh);
  status = (status == MPI_SUCCESS) ? 0 : 1;
#else /* HAVE_LIBLEMON */
  fh = fopen(filename, "r");
  status = (fh == NULL) ? 1 : 0;
  fflush(stderr);
#endif /* HAVE_LIBLEMON */

  if (status)
  {
    kill_with_error(fh, g_cart_id, "Could not open file. Aborting...\n");
  }

#ifdef HAVE_LIBLEMON
  *reader = lemonCreateReader(fh, g_cart_grid);
#else /* HAVE_LIBLEMON */
  *reader = limeCreateReader(fh);
#endif /* HAVE_LIBLEMON */

  if (*reader == (READER *)NULL)
    kill_with_error(fh, g_cart_id, "Could not create reader. Aborting...\n");
}

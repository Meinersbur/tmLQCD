
/* $Id$ */

/****************************************************
 * IO routines:
 *
 * write_gauge_field_time_p(char * filename)
 *   writes gauge field configuration to file
 *   with name filename.
 *
 * read_gauge_field_time_p(char * filename)
 *   reads gauge field configuration from file
 *   with name filename.
 *
 * Autor: Ines Wetzorke <wetzorke@ifh.de>
 *        Carsten Urbach <urbach@ifh.de>
 *
 ****************************************************/

#ifndef _IO_H
#define _IO_H

int write_gauge_field_time_p(char * filename);
int read_gauge_field_time_p(char * filename);
int write_spinorfield_eo_time_p(const int s, const int r, char * filename);
int read_spinorfield_eo_time(const int s, const int r, char * filename);

#endif

/***********************************************************************
* Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
*
* This file is part of tmLQCD.
*
* tmLQCD is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* tmLQCD is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#include "spinor.ih"

void write_propagator_format(LimeWriter *writer, paramsPropagatorFormat const *format) {
  uint64_t bytes;
  char *message;

  if(g_cart_id == 0) {
    message = (char*)malloc(512);
    sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	    "<etmcFormat>\n"
	    "  <field>diracFermion</field>\n"
	    "  <precision>%d</precision>\n"
	    "  <flavours>%d</flavours>\n"
	    "  <lx>%d</lx>\n"
	    "  <ly>%d</ly>\n"
	    "  <lz>%d</lz>\n"
	    "  <lt>%d</lt>\n"
	    "</etmcFormat>",
	    format->prec, format->flavours,
	    format->nx, format->ny, format->nx, format->nt);
    
    bytes = strlen(message);
    write_header(writer, 1, 1, "etmc-propagator-format", bytes);
    write_message(writer, message, bytes);
    limeWriterCloseRecord(writer);
    free(message);
  }
  return;
}



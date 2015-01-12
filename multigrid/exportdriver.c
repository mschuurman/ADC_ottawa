#include <dx/dx.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <signal.h>

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif
#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif

void driver(void) ;
#define MAX_STRLEN   200
#define MAX_CHANNELS 100

static int ordinal[MAX_CHANNELS] = { 0 } ;

int
main (int argc, char *argv[] )
{
  setlinebuf(stdout) ;
  setlinebuf(stderr) ;
  /*
   *  Main program (Fortran). It will come back to drvGridHead/drvGridData
   *  when some data is available for rendering.
   */
  driver() ;
  /*
   *  Shut down the connection to DX, and leave.
   */
  return 0 ;
  }

void
drvGridData(int *channel, int *points, float *corners, float *field, float *data, int *len_name, char *name)
{
	int	        npts ;
	Array  		values    = NULL ; /* Data points                   */
	Array  		coords    = NULL ; /* Regular grid                  */
	Array		mesh      = NULL ; /* Regular mesh                  */
        Array           efield    = NULL ; /* Array for electric field      */
        String          comment   = NULL ; /* Object comments               */
	Field  		total     = NULL ;
	float		delta [9] = { 0. } ;
	float		origin[3] ;
        double          datasize ;
	int		i, j ;
	char *		errmsg	= "" ;
        char            buf[MAX_STRLEN] ;
        char            filename[MAX_STRLEN] ;

  strncpy((void*)buf,name,min(*len_name,MAX_STRLEN)) ;
  buf[min(*len_name,MAX_STRLEN)] = 0 ;
  datasize = 2*sizeof(float)*((double)points[0])*((double)points[1])*((double)points[2]) ;
  if (datasize>=(2.*1024.*1024.*1024.-1024.)) {
    fprintf (stderr,"exportdriver: Passed %lf-byte object. OpenDX does not support objects exceeding 2GB\n", datasize) ;
    fprintf (stdout,"exportdriver: Passed %lf-byte object. OpenDX does not support objects exceeding 2GB\n", datasize) ;
    return ;
    }
  npts = points[0]*points[1]*points[2] ;
  /*
   *  Dropping out of this basic block means an error!
   */
  do {
    /*
     *  Create and initialize array for the data values on grid
     */
    errmsg = "create value array" ;
    values = DXNewArray(TYPE_FLOAT,CATEGORY_COMPLEX,0) ;
    if (!values) break ;
    errmsg = "fill value arrys with data - no memory?" ;
    if ( !DXAddArrayData(values,0,npts,(void *)data) ) break ;
    /* 
     * Set the dependency of the data to be on positions 
     */
    errmsg = "set dependency attribute" ;
    if (!DXSetStringAttribute((Object)values, "dep", "positions")) break ;
    /*
     *  Create object for the electric field direction and strength
     */
    errmsg = "create efield array" ;
    efield = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3) ;
    if (!efield) break ;
    errmsg = "fill efield with data - no memory?" ;
    if ( !DXAddArrayData(efield,0,1,field) ) break ;
    /*
     *  Build object comment
     */
    comment = DXNewString((void*)buf) ;
    /*
     *  Initialize position grid
     */
    for (i=0;i<3;i++) {
      origin[i]   =  corners[0+2*i] ;
      if (i==0) j = 0 ; 
      if (i==1) j = 4 ; 
      if (i==2) j = 8 ; 
      delta [j] = (corners[1+2*i] - corners[0+2*i])/max(1,points[i]-1) ;
      }
    errmsg = "create position grid" ;
    coords = DXMakeGridPositionsV(3,(void*)points,origin,delta) ;
    if (!coords) break ;
    /*
     *  Initialize connections mesh
     */
    errmsg = "create connection mesh" ;
    mesh   = DXMakeGridConnectionsV(3,(void*)points) ;
    if (!mesh) break ;
    /*
     *  Create field to hold the data
     */
    errmsg = "create empty field" ;
    total = DXNewField();
    if (!field) break ;
    /* 
     *  Drop out brand-new coordinates, values, and mesh into the 
     *  field. Zeroing out of the arrays after a successful assignment
     *  guarantees that the data won't be release twice in case
     *  of an error.
     */
    errmsg = "associate coordinates and values with a field" ;
    if ( !DXSetComponentValue(total, "data",        (Object)values ) ) break ;
    values = NULL ;
    if ( !DXSetComponentValue(total, "positions",   (Object)coords ) ) break ;
    coords = NULL ;
    if ( !DXSetComponentValue(total, "connections", (Object)mesh   ) ) break ;
    mesh   = NULL ;
    if ( !DXSetComponentValue(total, "comment",     (Object)comment) ) break ;
    comment= NULL ;
    if ( !DXSetComponentValue(total, "efield",      (Object)efield ) ) break ;
    efield = NULL ;
    /* 
     *  Finalize the field, and return it to DX. This includes 
     *  construction of the bounding box.
     */
    errmsg = "finalize field" ;
    if (!DXEndField(total)) break ;
    errmsg = "exporting field" ;
    snprintf(filename,MAX_STRLEN,"chan%02d-%06d.dx",*channel,ordinal[*channel]++) ;
    /*
     * Good for debugging
     * if (!DXExportDX((Object)total,filename,"dx text follows")) break ;
     */
    if (!DXExportDX((Object)total,filename,"dx")) break ;
    DXDelete((Object)total) ;
    return ;
    } while(0) ;
  /*
   * Error condition. Issue an error message, free memory, and return.
   */
  DXAddMessage(errmsg) ;
  fprintf (stderr,"exportdriver: %s\n",DXGetErrorMessage()) ;
  DXDelete((Object)total) ;
  DXDelete((Object)mesh) ;
  DXDelete((Object)values) ;
  DXDelete((Object)coords) ;
  DXDelete((Object)efield) ;
  DXDelete((Object)comment) ;
  }

void
drvGridHead(int *active)
{
  *active = 1 ;
  }

void
drvGridShow( void )
{
  }


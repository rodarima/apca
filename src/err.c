#include "err.h"
#include <stdio.h>
#include <stdlib.h>

/**  Error handler  **************************************************/

void erhand(char *err_msg)
{
    fprintf(stderr,"Run-time error:\n");
    fprintf(stderr,"%s\n", err_msg);
    fprintf(stderr,"Exiting to system.\n");
    exit(1);
}

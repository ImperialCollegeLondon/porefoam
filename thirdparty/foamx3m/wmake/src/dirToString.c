
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

int main()
{
    int c;
    int nextupper = 0;

    while ((c=getchar()) != EOF)
    {
        if (c == '/')
        {
            nextupper = 1;
        }
        else
        {
            if (nextupper)
            {
                putchar(toupper(c));
            }
            else
            {
                putchar(c);
            }

            nextupper = 0;
        }
    }

    return 0;
}


/*****************************************************************************/

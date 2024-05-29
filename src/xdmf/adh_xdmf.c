#include <stdio.h>
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Utitilities for XMF
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] g           (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void write_xdmf_header(FILE *xmf){

    fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
    //Domain information
    fprintf(xmf, "\t <Domain Name=\"Adh Sim\">\n");
    //Time based grid start
    fprintf(xmf, "\t\t<Grid Name=\"MeshTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
    

}

void write_xdmf_tail(FILE *xmf){

    fprintf(xmf, "\t\t\t</Grid>\n");
    //Grid Finish
    fprintf(xmf, "\t\t</Grid>\n");
    fprintf(xmf, "\t</Domain>\n");
    //Domain finished
    fprintf(xmf, "</Xdmf>\n");
    //XDMF File finished
    

}
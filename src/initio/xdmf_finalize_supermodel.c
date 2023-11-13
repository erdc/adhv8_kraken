
#include "global_header.h"
#ifdef _ADH_HDF5
#define HDF5_INIT
#include "adh_hdf5.h"
#define IN_LIBXML /* required for LibXML library configuration */
#include <libxml/tree.h>
#include <libxml/xpath.h>
#endif
/* ADH Version 2.0.0 6-04 */
/* Write master XMF file (if in parallel) and finish up XMDF/HDF5 */

#ifdef _ADH_HDF5
void xdmf_finalize_supermodel(SSUPER_MODEL *sm, int nsupermodels, int myid)
{
    if (myid > 0)  return;
    /* If we are serial, then we don't need to create the master XMF file */
    //if (npes <= 1) return;
    /* If we have just 1 model and 1 supermodel, then we don't need to create the master XMF file */
    if (nsupermodels == 1 && sm[0].nsubmodels == 1) return;
  
    int npes;
    SMODEL *mod=NULL;
    HDF5 *hdf5=NULL;
    herr_t status;
    xmlDtdPtr dtd = NULL;
    xmlNodePtr node1 = NULL;
    xmlNodePtr node2 = NULL;
    xmlNodePtr temporal_node = NULL;    /* the temporal COLLECTION node */
    xmlNodePtr spatial_node = NULL; /* the spatial COLLECTION node */
    xmlNodePtr root_node = NULL, root_node_surf=NULL;
    xmlDocPtr doc_out = NULL;     /* the concatenated document */
    xmlDocPtr doc_in = NULL;      /* pointer to individual processor document */
    xmlNodePtr current_collection = NULL; /* pointer in current doc to collection grid node */
    xmlNodePtr current_temporal = NULL;   /* pointer into current doc at level to copy */
    xmlXPathContextPtr xpathCtx;
    xmlChar *xpathExpr;
    xmlXPathObjectPtr xpathObj;
    xmlAttr *props = NULL;
    int i,j;
    int igrid = 0, isupmod=0, isubmod=0;
    int icount = 0;
    int num_time_steps = 0;
    char basename[128];
    char *proj_name;
    char filename[128];
    FILE *fp;
  
    xmlInitParser();
    sprintf(basename, "superxdmf");
   
    /*
     * Creates a new document, a node and set it as a root node
     */
    doc_out = xmlNewDoc(BAD_CAST "1.0");
    root_node = xmlNewNode(NULL, BAD_CAST "Xdmf");
    xmlNewProp(root_node, BAD_CAST "Version", BAD_CAST "2.0");
    xmlNewProp(root_node, BAD_CAST "xmlns:xi", BAD_CAST "http://www.w3.org/2001/XInclude");
    xmlDocSetRootElement(doc_out, root_node);
  
    /* loop over xmf files of all supermodels and submodels over all pes */
    for (isupmod=0; isupmod < nsupermodels; isupmod++){
        for (isubmod=0; isubmod < sm[isupmod].nsubmodels; isubmod++){
            mod = &(sm[isupmod].submodel[isubmod]);
            hdf5 = &(mod->hdf5);
            npes = mod->grid->smpi->npes;
            for (igrid = 0; igrid < npes; igrid++) {
                /* Advance to next processor's file */
                proj_name = mod->io->proj_name;
                sprintf(filename, "%s_p%d.xmf", proj_name, igrid);
                if (npes==1) sprintf(filename, "%s_pall1.xmf", proj_name);
                doc_in = xmlReadFile(filename, NULL, 0);
          
                /* Create xpath evaluation context */
                xpathCtx = xmlXPathNewContext(doc_in);
                if (xpathCtx == NULL) {
                    tl_error("XDMF Finalize Error: unable to create new XPath context\n");
                    xmlFree(doc_out);
                    exit(0);
                }
                xpathExpr = xmlCharStrdup("/Xdmf/Domain/Grid[@CollectionType=\"Temporal\"]");
                xpathObj = xmlXPathEvalExpression(xpathExpr, xpathCtx);
                if (xpathObj == NULL) {
                    tl_error("XDMF Finalize Error: unable to evaluate xpath expression ");
                    xmlXPathFreeContext(xpathCtx);
                    xmlFreeDoc(doc_out);
                    exit(0);
                }
          
                /* Gather all timestep nodes */
                current_collection = xpathObj->nodesetval->nodeTab[0];
          
                free(xpathExpr);
                free(xpathObj);
          
                /* the number of children at this point = no. of time steps */
                num_time_steps = xmlChildElementCount(current_collection);
          
                /* set first time step collection in tree */
                current_temporal = current_collection->children;
                while (current_temporal->type != XML_ELEMENT_NODE){
                    current_temporal = current_temporal->next;
                }
          
                /* Create file template */
                //if (igrid == 0) {
                if (igrid == 0 && isubmod == 0 && isupmod == 0) {
                    /* Creates a DTD declaration. Isn't mandatory. */
                    dtd = xmlCreateIntSubset(hdf5->doc, BAD_CAST "Xdmf", NULL, BAD_CAST "Xdmf.dtd");
          
                    /* xmlNewChild() creates a new node, which is "attached" as child node
                    * of root_node node. */
                    node1 = xmlNewChild(root_node, NULL, BAD_CAST "Domain", NULL);
          
                    /* GRID element for TEMPORAL collection */
                    temporal_node = xmlNewChild(node1, NULL, BAD_CAST "Grid", NULL);
                    xmlNewProp(temporal_node, BAD_CAST "CollectionType", BAD_CAST "Temporal");
                    xmlNewProp(temporal_node, BAD_CAST "GridType", BAD_CAST "Collection");
                    //xmlNewProp(temporal_node, BAD_CAST "Name", BAD_CAST "Mesh");
                    xmlNewProp(temporal_node, BAD_CAST "Name", BAD_CAST proj_name);
          
                    /* START OF SPATIAL COLLECTION */
                    for (i = 0; i < num_time_steps; i++) {
                        spatial_node = xmlNewChild(temporal_node, NULL, BAD_CAST "Grid", NULL);
                        xmlNewProp(spatial_node, BAD_CAST "CollectionType", BAD_CAST "Spatial");
                        xmlNewProp(spatial_node, BAD_CAST "GridType", BAD_CAST "Collection");
                    }
                }
          
                /* Start with first spatial node */
                spatial_node = temporal_node->children;
                while (spatial_node->type != XML_ELEMENT_NODE) {
                    spatial_node = spatial_node->next;
                }
          
                for (i = 0; i < num_time_steps; i++) {
                    /* Copy Time stamp above collection for XDMF animation in ParaView */
                    if (igrid == 0) {
                        node2 = xmlFirstElementChild(current_temporal);
                        node1 = xmlCopyNode(node2, 1);
                        xmlAddChild(spatial_node, node1);
                    }
          
                    node1 = xmlCopyNode(current_temporal, 1);
                    xmlAddChild(spatial_node, node1);
                    spatial_node = xmlNextElementSibling(spatial_node);
                    current_temporal = xmlNextElementSibling(current_temporal);
                }
                xmlFreeDoc(doc_in);
            }
        }
    }

    /* Create Master file */
    sprintf(filename,"%s_pall_%d.xmf", basename, igrid);
    xmlSaveFormatFileEnc(filename, doc_out, "UTF-8", 1);
  
    /*free the document */
    xmlFreeDoc(doc_out);
  
    /* loop over xmf files of all supermodels and submodels over all pes */
    int flag = 0;
    for (isupmod=0; isupmod < nsupermodels; isupmod++){
        for (isubmod=0; isubmod < sm[isupmod].nsubmodels; isubmod++){
            mod = &(sm[isupmod].submodel[isubmod]);
            if (mod->flag.SW3_FLOW){
                flag = 1;
                break;
            }
        }
    }
  
    /* Writing surface node data, only for SW3D models */

    if (flag==1) { /* mod->flag.SW3_FLOW */
        doc_out = xmlNewDoc(BAD_CAST "1.0");
        root_node_surf = xmlNewNode(NULL, BAD_CAST "Xdmf");
        xmlNewProp(root_node_surf, BAD_CAST "Version", BAD_CAST "2.0");
        xmlNewProp(root_node_surf, BAD_CAST "xmlns:xi", BAD_CAST "http://www.w3.org/2001/XInclude");
        xmlDocSetRootElement(doc_out, root_node_surf);
        for (isupmod=0; isupmod < nsupermodels; isupmod++){
            for (isubmod=0; isubmod < sm[isupmod].nsubmodels; isubmod++){
                mod = &(sm[isupmod].submodel[isubmod]);
                hdf5 = &(mod->hdf5);
                npes = mod->grid->smpi->npes;
                for (igrid = 0; igrid < npes; igrid++) {
                    /* Advance to next processor's file */
                    proj_name = mod->io->proj_name;
                    if (mod->flag.SW3_FLOW) {
                        sprintf(filename, "%s_surf_p%d.xmf", proj_name, igrid);
                        if (npes==1) sprintf(filename, "%s_surf_pall1.xmf", proj_name);
                    }
                    else if (mod->flag.SW2_FLOW){
                        sprintf(filename, "%s_p%d.xmf", proj_name, igrid);
                        if (npes==1) sprintf(filename, "%s_pall1.xmf", proj_name);
                    }
                    doc_in = xmlReadFile(filename, NULL, 0);
    
                    /* Create xpath evaluation context */
                    xpathCtx = xmlXPathNewContext(doc_in);
                    if (xpathCtx == NULL){
                        tl_error("XDMF Finalize Error: unable to create new XPath context\n");
                        xmlFree(doc_out);
                        exit(0);
                    }
                    xpathExpr = xmlCharStrdup("/Xdmf/Domain/Grid[@CollectionType=\"Temporal\"]");
                    xpathObj = xmlXPathEvalExpression(xpathExpr, xpathCtx);
                    if (xpathObj == NULL) {
                        tl_error("XDMF Finalize Error: unable to evaluate xpath expression ");
                        xmlXPathFreeContext(xpathCtx);
                        xmlFreeDoc(doc_out);
                        exit(0);
                    }
    
                    /* Gather all timestep nodes */
                    current_collection = xpathObj->nodesetval->nodeTab[0];
    
                    free(xpathExpr);
                    free(xpathObj);
    
                    /* the number of children at this point = no. of time steps */
                    num_time_steps = xmlChildElementCount(current_collection);
    
                    /* set first time step collection in tree */
                    current_temporal = current_collection->children;
                    while (current_temporal->type != XML_ELEMENT_NODE) {
                        current_temporal = current_temporal->next;
                    }
    
                    /* Create file template */
                    //if (igrid == 0) {
                    if (igrid == 0 && isubmod == 0 && isupmod == 0) {
                        /* Creates a DTD declaration. Isn't mandatory. */
                        dtd = xmlCreateIntSubset(hdf5->doc_surf, BAD_CAST "Xdmf", NULL, BAD_CAST "Xdmf.dtd");
    
                        /* xmlNewChild() creates a new node, which is "attached" as child node
                        * of root_node node. */
                        node1 = xmlNewChild(root_node_surf, NULL, BAD_CAST "Domain", NULL);
    
                        /* GRID element for TEMPORAL collection */
                        temporal_node = xmlNewChild(node1, NULL, BAD_CAST "Grid", NULL);
                        xmlNewProp(temporal_node, BAD_CAST "CollectionType", BAD_CAST "Temporal");
                        xmlNewProp(temporal_node, BAD_CAST "GridType", BAD_CAST "Collection");
                        xmlNewProp(temporal_node, BAD_CAST "Name", BAD_CAST proj_name);
    
                        /* START OF SPATIAL COLLECTION */
                        for (i = 0; i < num_time_steps; i++) {
                            spatial_node = xmlNewChild(temporal_node, NULL, BAD_CAST "Grid", NULL);
                            xmlNewProp(spatial_node, BAD_CAST "CollectionType", BAD_CAST "Spatial");
                            xmlNewProp(spatial_node, BAD_CAST "GridType", BAD_CAST "Collection");
                        }
                    }
    
                    /* Start with first spatial node */
                    spatial_node = temporal_node->children;
                    while (spatial_node->type != XML_ELEMENT_NODE) {
                        spatial_node = spatial_node->next;
                    }
    
                    for (i = 0; i < num_time_steps; i++) {
                        /* Copy Time stamp above collection for XDMF animation in ParaView */
                        if (igrid == 0) {
                            node2 = xmlFirstElementChild(current_temporal);
                            node1 = xmlCopyNode(node2, 1);
                            xmlAddChild(spatial_node, node1);
                        }
    
                        node1 = xmlCopyNode(current_temporal, 1);
                        xmlAddChild(spatial_node, node1);
                        spatial_node = xmlNextElementSibling(spatial_node);
                        current_temporal = xmlNextElementSibling(current_temporal);
                    }
                    xmlFreeDoc(doc_in);
                }
            }
        }
        /* Create Master file */
        sprintf(filename,"%s_surf_pall_%d.xmf", basename, igrid);
        xmlSaveFormatFileEnc(filename, doc_out, "UTF-8", 1);
    
        /*free the document */
        xmlFreeDoc(doc_out);
    }
  
    return;
}
#endif


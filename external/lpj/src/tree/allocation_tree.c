/***************************************************************************/
/**                                                                       **/
/**        a  l  l  o  c  a  t  i  o  n  _  t  r  e  e  .  c              **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 07.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"

#define CDEBT_MAXLOAN_DEFICIT 0.8 /* maximum loan as a fraction of deficit*/
#define CDEBT_MAXLOAN_MASS 0.2 /* maximum loan as a fraction of (sapwood-cdebt)*/
#define NSEG 20	/* number of segments (parameter in numerical methods)*/

typedef struct
{
  Real k1,lm,k3,b,ind_leaf,ind_heart;
} Data;

static Real f(Real leaf_inc,Data *data) 
{
 return  data->k1*(data->b-leaf_inc*data->lm+data->ind_heart)-
	 pow((data->b-leaf_inc*data->lm)/(data->ind_leaf+leaf_inc)*data->k3,
             1.0+2/allom3);
} /* of 'f' */

void allocation_tree(Litter *litter,Pft *pft,Real *fpc_inc)
{
  Real bm_inc_ind,lmtorm;
  Real cmass_loan,cmass_deficit;
  Real x1,x2,swap;
  Treephys2 tinc_ind,tinc_ind_min;
  Pfttree *tree;
  Data data; 
  bm_inc_ind=pft->bm_inc/pft->nind;
  lmtorm=getpftpar(pft,lmro_ratio)*(pft->wscal_mean/365);
  tree=pft->data;
  tinc_ind.heartwood=tinc_ind.debt=0.0;
  
  if (lmtorm<1.0e-10) 	/* No leaf production possible - put all biomass 
                           into roots (Individual will die next time period)*/
  {
    tinc_ind.leaf=0.0;
    tinc_ind.root=bm_inc_ind;
    tinc_ind.sapwood=-tree->ind.sapwood;
    tinc_ind.heartwood=-tinc_ind.sapwood;
    tinc_ind.debt=0.0;
  }
  else
  { 
    if (tree->height>0.0) 
    {
      tinc_ind_min.leaf=k_latosa*tree->ind.sapwood/(wooddens*tree->height*
                        pft->par->sla)-tree->ind.leaf;
      tinc_ind_min.root=k_latosa*tree->ind.sapwood/(wooddens*tree->height*
                        pft->par->sla*lmtorm)-tree->ind.root;
    }
    else 
      tinc_ind_min.leaf=tinc_ind_min.root=0.0;

    cmass_deficit=tinc_ind_min.leaf+tinc_ind_min.root-bm_inc_ind;
    if (cmass_deficit>0.0) 
    {
      cmass_loan=max(min(cmass_deficit*CDEBT_MAXLOAN_DEFICIT,
                 tree->ind.sapwood-tree->ind.debt)*CDEBT_MAXLOAN_MASS,0.0);
      bm_inc_ind+=cmass_loan;
      tinc_ind.debt=cmass_loan;
    }
    else 
      tinc_ind.debt=0.0;
  
    if (tinc_ind_min.root>=0.0 && tinc_ind_min.leaf>=0.0 &&
        tinc_ind_min.root+tinc_ind_min.leaf<=bm_inc_ind || bm_inc_ind<=0.0) 
    {
      data.b= tree->ind.sapwood+bm_inc_ind-tree->ind.leaf/lmtorm+
              tree->ind.root;
      data.lm=1+1/lmtorm;
      data.k1=pow(allom2,2.0/allom3)*4.0*M_1_PI/wooddens;
      data.k3=k_latosa/wooddens/pft->par->sla;
      data.ind_leaf=tree->ind.leaf;
      data.ind_heart=tree->ind.heartwood;
      x2=(bm_inc_ind-(tree->ind.leaf/lmtorm-tree->ind.root))/data.lm;
      x1= (tree->ind.leaf<1.0e-10)  ? x2/NSEG : 0;

/*		Bisection loop
 *		Search iterates on value of xmid until xmid lies within
 *		xacc of the root, i.e. until |xmid-x|<xacc where f(x)=0
 */
 
      if(x1==0 && x2==0 || data.b-x1*data.lm<0.0 || data.ind_leaf+x1<=0.0 
         || data.b-x2*data.lm<0.0 || data.ind_leaf+x2<=0.0 )
        tinc_ind.leaf=0;
      else
        tinc_ind.leaf=leftmostzero((Bisectfcn)f,x1,x2,&data,0.001,1.0e-10,40);
      if (tinc_ind.leaf<0.0)
       tinc_ind.root=0.0;
      else
       tinc_ind.root=(tinc_ind.leaf+tree->ind.leaf)/lmtorm-tree->ind.root; 
      tinc_ind.sapwood=bm_inc_ind-tinc_ind.leaf-tinc_ind.root;
    }
    else 
    {
    
/* 		Abnormal allocation: reduction in some biomass compartment(s) to
*  		satisfy allometry
*  		Attempt to distribute this year's production among leaves and roots only
*/
      tinc_ind.leaf=(bm_inc_ind+tree->ind.root-tree->ind.leaf/lmtorm)/
                    (1.0+1.0/lmtorm);
      if (tinc_ind.leaf>0.0) 
      {
        tinc_ind.root=bm_inc_ind-tinc_ind.leaf;
      }
      else 
      {
        tinc_ind.root=bm_inc_ind;
        tinc_ind.leaf=(tree->ind.root+tinc_ind.root)*lmtorm-tree->ind.leaf;
        litter->ag_tree+=-tinc_ind.leaf*pft->nind;
      }
      tinc_ind.sapwood=(tinc_ind.leaf+tree->ind.leaf)*wooddens*tree->height*
                       pft->par->sla/k_latosa-tree->ind.sapwood;
      tinc_ind.heartwood=-tinc_ind.sapwood;
    }
  }
  tree->ind.leaf+=tinc_ind.leaf;
  tree->ind.root+=tinc_ind.root;
  tree->ind.sapwood+=tinc_ind.sapwood;
  tree->ind.heartwood+=tinc_ind.heartwood;
  tree->ind.debt+=tinc_ind.debt;
#ifdef DEBUG
  printf ("leaf %f root %f sapwood %f heartwood %f\n",
    tinc_ind.leaf,tinc_ind.root,tinc_ind.sapwood,tinc_ind.heartwood);
  printf ("debt %f debt_ind %f sapwood_ind %f bm_inc %f bm_inc_ind %f\n",
    tinc_ind.debt,tree->ind.debt,tree->ind.sapwood,pft->bm_inc/pft->nind,bm_inc_ind);
#endif
  allometry_tree(pft);
  *fpc_inc=fpc_tree(pft);
} /* of 'allocation_tree' */

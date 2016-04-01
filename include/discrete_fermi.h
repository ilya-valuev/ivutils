/*e***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2014        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD
 *
 *   $Revision: 1.2 $
 *   $Date: 2014/10/28 10:50:48 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/discrete_fermi.h,v 1.2 2014/10/28 10:50:48 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/discrete_fermi.h,v $
$Revision: 1.2 $
$Author: valuev $
$Date: 2014/10/28 10:50:48 $
*/
/*e****************************************************************************
 * $Log: discrete_fermi.h,v $
 * Revision 1.2  2014/10/28 10:50:48  valuev
 * adapting tcpengine for AWPMD
 *
 * Revision 1.1  2014/07/28 08:39:49  valuev
 * discrete fermi distributions
 *
 *
*******************************************************************************/
#ifndef __DISCRETE_FERMI_H
#define __DISCRETE_FERMI_H

/**\en @file discrete_fermi.h
       @brief Functions for calculating discrete Fermi distributions in multi-level multi electron systems.
        
**/

# include <stdio.h>

double exc_en_fermi(int nlev, int ne, double beta_de, FILE *f=NULL, double equant = .1); 

double exc_en_hartree(int nlev, int ne, double beta_de, FILE *f=NULL, double equant = .1);

/*

double exc_en_discr_fermi(const vector<double> &levels,const vector<int> &g, size_t ne, double beta, vector<double> &occupancy){
  vector<double> factors(levels.size());
  for(size_t i=0;i<levels.size();i++)
    factors[i] = exp(-levels[i]*beta);
}*/

#endif

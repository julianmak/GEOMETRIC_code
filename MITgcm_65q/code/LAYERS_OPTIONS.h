C $Header: /u/gcmpack/MITgcm/pkg/layers/LAYERS_OPTIONS.h,v 1.1 2009/09/16 21:25:47 rpa Exp $
C $Name: checkpoint62x $

C
C Use this file for selecting options within package "LAYERS"

#ifndef LAYERS_OPTIONS_H
#define LAYERS_OPTIONS_H

#include "PACKAGES_CONFIG.h"

#ifdef ALLOW_LAYERS

#include "CPP_OPTIONS.h"

C Compute isopycnal tranports in the U direction?
#define LAYERS_UFLUX
C Compute isopycnal tranports in the V direction?
#define LAYERS_VFLUX
C Keep track of layer thicknesses and the velocity in each layer?
#define LAYERS_THICKNESS

C The MNC stuff is too complicated
#undef LAYERS_MNC

C Allow use of potential density as a layering field.
#define LAYERS_PRHO_REF

#endif /* ALLOW_LAYERS */
#endif /* LAYERS_OPTIONS_H */


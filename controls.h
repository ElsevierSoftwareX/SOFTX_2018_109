#ifndef CONTROLS_H
#define CONTROLS_H

#include "datatypes.h"

typedef struct MainControls_ MainControls;
typedef struct ContactControls_ ContactControls;
typedef struct AreaControls_ AreaControls;

typedef struct OPControls_ OPControls;
typedef struct OPUncControls_ OPUncControls;
typedef struct OPMCControls_ OPMCControls;

typedef struct OPODRControls_ OPODRControls;
typedef struct OPODRUncControls_ OPODRUncControls;
typedef struct OPODRMCControls_ OPODRMCControls;
typedef struct OPODRUnc2Controls_ OPODRUnc2Controls;

typedef struct TangentControls_ TangentControls;
typedef struct TangentUncControls_ TangentUncControls;
typedef struct TangentMCControls_ TangentMCControls;

typedef struct HertzControls_ HertzControls;
typedef struct HertzUncControls_ HertzUncControls;
typedef struct HertzMCControls_ HertzMCControls;

typedef struct HertzODRControls_ HertzODRControls;
typedef struct HertzODRUncControls_ HertzODRUncControls;
typedef struct HertzODRMCControls_ HertzODRMCControls;

typedef struct StiffnessControls_ StiffnessControls;
typedef struct StiffnessUncControls_ StiffnessUncControls;
typedef struct StiffnessMCControls_ StiffnessMCControls;

typedef struct SlopesControls_ SlopesControls;
typedef struct SlopesUncControls_ SlopesUncControls;
typedef struct SlopesMCControls_ SlopesMCControls;

typedef struct Ph2Controls_ Ph2Controls;
typedef struct ApopinControls_ ApopinControls;

typedef struct WorkControls_ WorkControls;
typedef struct WorkUncControls_ WorkUncControls;


typedef struct {
    MainControls *maincontrols;
    ContactControls *contactcontrols;
    AreaControls *areacontrols;

    OPControls *opcontrols;
    OPUncControls *opunccontrols;
    OPMCControls *opmccontrols;

    OPODRControls *opodrcontrols;
    OPODRUncControls *opodrunccontrols;
    OPODRMCControls *opodrmccontrols;

    TangentControls *tgcontrols;
    TangentUncControls *tgunccontrols;
    TangentMCControls *tgmccontrols;

    HertzControls *hertzcontrols;
    HertzUncControls *hertzunccontrols;
    HertzMCControls *hertzmccontrols;

    HertzODRControls *hertzodrcontrols;
    HertzODRUncControls *hertzodrunccontrols;
    HertzODRMCControls *hertzodrmccontrols;

    StiffnessControls *stiffnesscontrols;
    StiffnessUncControls *stiffnessunccontrols;
    StiffnessMCControls *stiffnessmccontrols;

    SlopesControls *slopescontrols;
    SlopesUncControls *slopesunccontrols;
    SlopesMCControls *slopesmccontrols;

    Ph2Controls *ph2controls;

    ApopinControls *apopcontrols;

    WorkControls *workcontrols;
    WorkUncControls *workunccontrols;
} Controls;

typedef struct {
    Args *args;
    Controls *controls;
    struct Functions *functions;
} Data;

#endif

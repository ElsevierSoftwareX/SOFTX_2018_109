#ifndef DATATYPES_H
#define DATATYPES_H


/* Area */
#include "instrument.h"

#include <libprocess/gwyprocess.h>

/* TEMPORARY SOLUTION UNTIL TYPE DEFINITIONS ARE DISTRIBUTED TO THE TOOLS */

typedef struct {
    gdouble slope;
    gdouble intercept;
    gdouble chi; //kdo vi co toje
    gdouble corr; //kdo vi co toje
    gdouble R2, R2adj, chi2;
    gdouble SSres;
    gdouble unc_slope;
    gdouble unc_intercept;
    gdouble cov_slope_intercept;
    gdouble unc_slope_x;
    gdouble unc_slope_y;
    gdouble unc_intercept_x;
    gdouble unc_intercept_y;
    gdouble cov_slope_intercept_x;
    gdouble cov_slope_intercept_y;
    gdouble *auxtoberemoved_axi; // TODO remove this asap
} LinReg;


enum PointSelStatus {
    POINT_STATUS_NOT_SET = 0,
    POINT_STATUS_SET_AUTO,
    POINT_STATUS_SET_MAN,
    POINT_STATUS_FROM_FILE
};


/* pasted back from controls.h */
#define ndatafileformats 7

enum DataFileFormat {
    DATA_FILE_FORMAT_NIGET = 0,
    DATA_FILE_FORMAT_UNHT_D_NM_F_MN,
    DATA_FILE_FORMAT_UNHT_D_NM_F_UN,
    DATA_FILE_FORMAT_HYSITRON,
    DATA_FILE_FORMAT_D_NM_F_MN,
    DATA_FILE_FORMAT_F_MN_D_NM,
    DATA_FILE_FORMAT_F_UN_D_NM
};
/* */

enum ExportFormat {EXPORT_FORMAT_PLAIN = 0
                  };

typedef struct {
    gchar *filename;
    enum DataFileFormat fileformat;

    /* all h in nm, all F in mN */
    GwyDataLine *horig;
    GwyDataLine *Forig;
    GwyDataLine *hload;
    GwyDataLine *Fload;
    GwyDataLine *hhold;
    GwyDataLine *Fhold;
    GwyDataLine *hunload;
    GwyDataLine *Funload;

    GwyDataLine *time;
    GwyDataLine *htime;
    GwyDataLine *Ftime;

    gint ndata;
    gdouble hmax;
    gdouble Fmax;
    gdouble hc;
    gdouble Fc;

    /* with respect to horig, Forig */
    gint i_contact_load;
    gint i_load_hold;
    gint i_hold_unload;
    gint i_end_unload;

    enum PointSelStatus status_contact_load, status_load_hold,
         status_hold_unload, status_end_unload;
} FDdata;

/* OP */

typedef struct {
    gdouble hp, hc, hr, eps, m, S, Aphc, alpha; /* hp, hc, hr in nm, S in mN/nm, Aphc in nm^2 */
    gdouble Eit, Er, Hit; /* Eit, Er in Pa, Hit in MPa */
    gdouble beta;
    LinReg reghp;
    LinReg regS;
    gboolean has_fit_hp, has_fit_S;
    gdouble hpfrom, hpto;
    gdouble hpfrom_pct_Fmax, hpto_pct_Fmax;
    gdouble hprange[2];
    gdouble hprange_pct_Fmax[2];
    gdouble Sfrom, Sto;
    gdouble Srange[2];
    gdouble Srange_pct_Fmax[2];
    gdouble Sfrom_pct_Fmax, Sto_pct_Fmax;
    GwyDataLine *xhpfit, *yhpfit;
    GwyDataLine *xSfit, *ySfit;
    gboolean Finputhp, FinputS;
    gint nfithpdata, nfitSdata;
    //hpfrom,hpto, Sfrom,Sto  are the data from the entries, i.e. from the selection
    //hprange, Srange are the ranges that were actually used in the fit
} OPdata;

typedef struct {
    gdouble uh, uF;
    /*    gdouble unu, unui, uEi; */

    gint Nmc; // Monte Carlo iterations

    gdouble uSuh, uSuF;
    gdouble uhcuh, uhcuF;
    gdouble uAuh, uAuF;
    gdouble uHituh, uHituF;

    gdouble uEruh, uEruF;
    gdouble uEituh, uEituF;
    gdouble uEitunu, uEitunui, uEituEi;

    gdouble uS, uhc, uA, uHit, uEr, uEit;

    gdouble *w;
    gdouble *Sc, *Ec, *Hc;
} OPUncdata;

typedef struct {
    OPdata *opmc;
    gint Nmc; // Monte Carlo iterations
    gint skipmc; //how many iterations in Monte Carlo are nonsense
    gdouble *mcavg, *mcstd;
    gdouble **mcxhist, * *mcyhist;
    gint mcnstat;
    gdouble **mcdd;
} OPMCdata;

/* OP ODR */

typedef struct {
    gdouble hp, hc, hr, eps, m, S, Aphc, alpha; /* hp, hc, hr in nm, S in mN/nm, Aphc in nm^2 */
    gdouble Eit, Er, Hit; /* Eit, Er in Pa, Hit in MPa */
    gdouble SSres, R2, R2adj, chi2;
    gdouble beta;
    gboolean has_fit;
    gdouble from, to;
    gdouble from_pct_Fmax, to_pct_Fmax;
    gdouble range[2];
    gdouble range_pct_Fmax[2];
    gboolean radial_corr;
    gdouble radial_angle;
    GwyDataLine *xfit, *yfit;
    gboolean Finput;
    gint nfitdata;
    gchar *logfnm;
    gint infolog;
    gint fitinfo;
} OPODRdata;

typedef struct {

    gint Nmc; // Monte Carlo iterations

    gdouble uhp, covmhp;
    gdouble ueps;

    // contributions from uh,uF
    gdouble uhcuh, uhcuF;
    gdouble uAuh, uAuF;
    gdouble uSuh, uSuF;
    gdouble umuh, umuF;
    gdouble uHituh, uHituF;
    gdouble uEruh, uEruF;
    gdouble uEituh, uEituF;

    //noise contributions
    gdouble uAnoise;
    gdouble uHitnoise;
    gdouble uErnoise;
    gdouble uEitnoise;

    //area coeff. contributions
    gdouble uAcoeff;
    gdouble uHitcoeff;
    gdouble uErcoeff;
    gdouble uEitcoeff;

    //radial contributions
    gdouble uAradial;
    gdouble uHitradial;
    gdouble uErradial;
    gdouble uEitradial;

    //sample-indenter contributions
    gdouble uEitunu, uEitunui, uEituEi;

    //total contributions
    gdouble uS, um, uhc,  uHit, uEr, uEit;
    gdouble uA;

    //shift of contact
    gdouble *w;
    gdouble *Ec, *Hc, *Ac, *hc;

    Instdata instdata;
} OPODRUncdata;

typedef struct {
    OPODRdata *opodrmc;
    gint Nmc; // Monte Carlo iterations
    gint skipmc; //how many iterations in Monte Carlo are nonsense
    gdouble *mcavg, *mcstd;
    gdouble **mcxhist, * *mcyhist;
    gint mcnstat;
    gdouble **mcdd;
} OPODRMCdata;

/* Tangent */

typedef struct {
    gdouble hc, S, Aphc;/* hp in nm, S in mN/nm, Aphc in nm^2 */
    gdouble eps; // =0.75
    gdouble beta;
    LinReg reg;
    gdouble Eit, Er, Hit; /* Eit, Er in Pa, Hit in MPa */
    gboolean has_fit;
    gdouble from, to;
    gdouble range[2];
    gdouble from_pct_Fmax, to_pct_Fmax;
    gdouble range_pct_Fmax[2];
    GwyDataLine *xfit, *yfit;
    gboolean Finput;
    gint nfitdata;
} Tangentdata;

typedef struct {
    gdouble uh, uF;
    /*    gdouble unu, unui, uEi; */

    gdouble uSuh, uSuF;
    gdouble uhcuh, uhcuF;
    gdouble uAuh, uAuF;
    gdouble uHituh, uHituF;

    gdouble uEruh, uEruF;
    gdouble uEituh, uEituF;
    gdouble uEitunu, uEitunui, uEituEi;

    gdouble uS, uhc, uA, uHit, uEr, uEit;

    gdouble *w;
    gdouble *Sc, *Ec, *Hc;

    gint Nmc;
} TangentUncdata;

typedef struct {
    Tangentdata *tgmc;
    gint Nmc;
    gint skipmc; //how many iterations in Monte Carlo are nonsense
    gdouble *mcavg, *mcstd;
    gdouble **mcxhist, * *mcyhist;
    gint mcnstat;
    gdouble **mcdd;
} TangentMCdata;

/* Hertz */

typedef struct {
    gdouble radius; /* in m */
    LinReg reg;
    gint mode; /* 0 input R, 1 input Er, 2 input Eit */
    gdouble Eit, Er; /* in GPa */
    gboolean has_fit;
    gdouble from, to;
    gdouble range[2];
    gint nfitdata;
    GwyDataLine *xfit, *yfit;
} Hertzdata;


typedef struct {
    gdouble uh, uF;
    /*    gdouble unu, unui, uEi; */

    /* Rmode uradius input, uErtotal and uEittotal output
     * Ermode uEr input, uradiustotal output
     * Eitmode uEr input, uradiustotal output */
    gdouble uradius;
    gdouble uEr, uEit;

    /* Rmode */
    gdouble uEruh, uEruF;
    gdouble uEruradius;
    gdouble uErtotal;

    gdouble uEituh, uEituF;
    gdouble uEituradius;
    gdouble uEitunu, uEitunui, uEituEi;
    gdouble uEittotal;

    /* Ermode, Eitmode */
    gdouble uradiusuh, uradiusuF;
    gdouble uradiusuEr, uradiusuEit;
    gdouble uradiusunu, uradiusunui, uradiusuEi;
    gdouble uradiustotal;

    gdouble *w;
    gdouble *ERc; /*either Er or radius, depending on mode */

    gint Nmc; // Monte Carlo iterations
} HertzUncdata;

typedef struct {
    Hertzdata *hertzmc;
    gint Nmc;
    gint skipmc; //how many iterations in Monte Carlo are nonsense
    gdouble *mcavg, *mcstd;
    gdouble **mcxhist, * *mcyhist;
    gint mcnstat;
    gdouble **mcdd;
} HertzMCdata;

/* Hertz ODR */
typedef struct {
    gdouble radius; /* in m */
    gint mode; /* 0 input R, 1 input Er, 2 input Eit */
    gdouble Eit, Er; /* in GPa */
    gdouble h0, n, gamma; /*fit by function gamma(h-h0)^n */
    gint ifixb[3]; /*which parameters should be fixed: 0 - fixed, 1 - fitted */
    gdouble SSres, R2, R2adj, chi2; /*quality of fit */
    gboolean has_fit;
    gdouble from, to;
    gdouble range[2];
    gint nfitdata;
    GwyDataLine *xfit, *yfit;
    gchar *logfnm;
    gint infolog;
    gint fitinfo;
    gboolean radial_corr;
} HertzODRdata;

typedef struct {

    gint Nmc; // Monte Carlo iterations
    gdouble uh0uh, uh0uF;
    gdouble unuh, unuF;
    gdouble ugammauh, ugammauF;
    gdouble ugamma, un, uh0;

    /* Rmode uradius input, uErtotal and uEittotal output
         * Ermode uEr input, uradiustotal output
         * Eitmode uEr input, uradiustotal output */
    gdouble uradius;
    gdouble uEr, uEit;

    /* Rmode */
    gdouble uEruh, uEruF;
    gdouble uEruradius;
    gdouble uErtotal;

    gdouble uEituh, uEituF;
    gdouble uEituradius;
    gdouble uEitunu, uEitunui, uEituEi;
    gdouble uEittotal;

    /* Ermode, Eitmode */
    gdouble uradiusuh, uradiusuF;
    gdouble uradiusuEr, uradiusuEit;
    gdouble uradiusunu, uradiusunui, uradiusuEi;
    gdouble uradiustotal;

    gdouble *w;
    gdouble *ERc; /*either Er or radius, depending on mode */

    Instdata instdata;
} HertzODRUncdata;

typedef struct {
    HertzODRdata *hertzodrmc;
    gint Nmc;
    gint skipmc; //how many iterations in Monte Carlo are nonsense
    gdouble *mcavg, *mcstd;
    gdouble **mcxhist, * *mcyhist;
    gint mcnstat;
    gdouble **mcdd;
} HertzODRMCdata;

/* Stiffness */

typedef struct {
    /*loading fit  F = kl*h + ql */
    gdouble kl, ql;
    /*unloading fit  F = ku*h + qu */
    gdouble ku, qu;


    gboolean has_fit_load, has_fit_unload;

    /*ranges */
    //loadfrom,loadto, unloadfrom,unloadto  are the data from the entries, i.e. from the selection
    //loadrange, unloadrange are the ranges that were actually used in the fit
    gdouble loadfrom, loadto;
    gdouble loadrange[2];
    gdouble loadSSres, loadR2, loadR2adj, loadchi2; /*quality of fit */

    gdouble unloadfrom, unloadto;
    gdouble unloadrange[2];
    gdouble unloadSSres, unloadR2, unloadR2adj, unloadchi2; /*quality of fit */

    GwyDataLine *xloadfit, *yloadfit;
    GwyDataLine *xunloadfit, *yunloadfit;
    gint nfitloaddata, nfitunloaddata;

    gchar *logfnmload;
    gchar *logfnmunload;
    gint infologload;
    gint infologunload;
    gint fitinfoload;
    gint fitinfounload;
} Stiffnessdata;

typedef struct {

    gint Nmc; // Monte Carlo iterations

    gdouble ukl, uql;
    gdouble uku, uqu;
    gdouble ukluh, ukluF;
    gdouble ukuuh, ukuuF;
    gdouble uqluh, uqluF;
    gdouble uquuh, uquuF;

    gdouble *w;
    gdouble *klc, *qlc, *kuc, *quc;

    Instdata instdata;
} StiffnessUncdata;

typedef struct {
    Stiffnessdata *stiffnessmc;
    gint Nmc; // Monte Carlo iterations
    gint skipmc; //how many iterations in Monte Carlo are nonsense
    gdouble *mcavg, *mcstd;
    gdouble **mcxhist, * *mcyhist;
    gint mcnstat;
    gdouble **mcdd;
} StiffnessMCdata;



/* Two slopes */

typedef struct {
    /*loading fit  F = gamma (h-h0) ^n */
    gdouble gamma, h0, n; /* h0 in nm, n number, gamma weirs */
    /* unloading fit  F = alpha (h-hp) ^m */
    gdouble alpha, hp, m; /* hp in nm, m number, alpha weirs */

    gdouble eps, Sload, Sunload, Aphc; /*S in mN/nm, Aphc in nm^2, eps number */
    gdouble beta;
    gdouble Eit, Er, Hit; /* Eit, Er in Pa, Hit in MPa */

    /* LinReg regload;
        LinReg regunload;  */
    gboolean has_fit_load, has_fit_unload;

    /*ranges */
    //loadfrom,loadto, unloadfrom,unloadto  are the data from the entries, i.e. from the selection
    //loadprange, unloadrange are the ranges that were actually used in the fit
    gdouble loadfrom, loadto;
    gdouble loadfrom_pct_Fmax, loadto_pct_Fmax;
    gdouble loadrange[2];
    gdouble loadrange_pct_Fmax[2];
    gint loadifixb[3]; /*which parameters should be fixed: 0 - fixed, 1 - fitted */
    gdouble loadSSres, loadR2, loadR2adj, loadchi2; /*quality of fit */

    gdouble unloadfrom, unloadto;
    gdouble unloadfrom_pct_Fmax, unloadto_pct_Fmax;
    gdouble unloadrange[2];
    gdouble unloadrange_pct_Fmax[2];
    gint unloadifixb[3]; /*which parameters should be fixed: 0 - fixed, 1 - fitted */
    gdouble unloadSSres, unloadR2, unloadR2adj, unloadchi2; /*quality of fit */

    GwyDataLine *xloadfit, *yloadfit;
    GwyDataLine *xunloadfit, *yunloadfit;
    gboolean Finputload, Finputunload;
    gint nfitloaddata, nfitunloaddata;

    gchar *logfnmload;
    gchar *logfnmunload;
    gint infologload;
    gint infologunload;
    gint fitinfoload;
    gint fitinfounload;
} Slopesdata;

typedef struct {
    /*    gdouble unu, unui, uEi; */

    gint Nmc; // Monte Carlo iterations

    gdouble uSluh, uSluF;
    gdouble uSuuh, uSuuF;
    gdouble umuh, umuF;
    gdouble unuh, unuF;

    gdouble uAuh, uAuF;
    gdouble uHituh, uHituF;

    gdouble uEruh, uEruF;
    gdouble uEituh, uEituF;
    gdouble uEitunu, uEitunui, uEituEi;

    gdouble uSl, uSu, um, un, uA, uHit, uEr, uEit;

    gdouble *w;
    gdouble *Ec, *Hc, *Ac;

    Instdata instdata;
} SlopesUncdata;

typedef struct {
    Slopesdata *slopesmc;
    gint Nmc; // Monte Carlo iterations
    gint skipmc; //how many iterations in Monte Carlo are nonsense
    gdouble *mcavg, *mcstd;
    gdouble **mcxhist, * *mcyhist;
    gint mcnstat;
    gdouble **mcdd;
} SlopesMCdata;


/* Ph2 */

typedef struct {
    gint nmove;
    GwyDataLine *havg; /* in nm */
    GwyDataLine *Favg; /* in mN */
    GwyDataLine *Ph2;  /* in MPa */
    GwyDataLine *dPdh2;  /* in MPa */
    GwyDataLine *ind_remove;
} Ph2data;

/* Pop-ins A */

typedef struct {
    gint nmove;
    gint npopin;
    gdouble thresh, thresh2; // threshold for derivative to identify popin and to determine its beginning and end
    gdouble hpop; //minimum popin height
    gint wpop; // minimum popin width (number of datapoints in popin)
    gdouble *Fpopin; // Fpopin[npopin]
    gdouble *hpopin; // hpopin[npopin]
    gdouble *dhpopin;// dhpopin[npopin]
    gint *ileft; // index of left end of popin
    gint *iright; // index of right end of popin
    GwyDataLine *havg; // sliding average of lead
    GwyDataLine *dh; // derivative dh/di
    GwyDataLine *t;  // index i
    GwyDataLine *dth; // constant thresh
    GwyDataLine *dth2; //constant thresh2
} Apopindata;

/* Work */

typedef struct {
    /* all work in pJ */
    gdouble worktotal;
    gdouble workelastic;
    gdouble workplastic;
    gdouble eta;  // = W_e/W_t = W_e/(W_p + W_e)
    gint nmove; // width of moving average window
    GwyDataLine *hloadavg, *Floadavg; //moving averages
    GwyDataLine *hholdavg, *Fholdavg; //moving averages
    GwyDataLine *hunloadavg, *Funloadavg; //moving averages
} Workdata;

typedef struct {
    gdouble *We, *Wp;
} WorkUncdata;

/* Args */

typedef struct {
    FDdata  fddata;
    Area area;
    Instdata instdata;

    OPdata opdata;
    OPUncdata opunc;
    OPMCdata opmc;

    OPODRdata opodrdata;
    OPODRUncdata opodrunc;
    OPODRMCdata opodrmc;

    Tangentdata tgdata;
    TangentUncdata tgunc;
    TangentMCdata tgmc;

    Hertzdata hertzdata;
    HertzUncdata hertzunc;
    HertzMCdata hertzmc;

    HertzODRdata hertzodrdata;
    HertzODRUncdata hertzodrunc;
    HertzODRMCdata hertzodrmc;

    Stiffnessdata stiffnessdata;
    StiffnessUncdata stiffnessunc;
    StiffnessMCdata stiffnessmc;

    Slopesdata slopesdata;
    SlopesUncdata slopesunc;
    SlopesMCdata slopesmc;

    Ph2data ph2data;
    Apopindata apopdata;

    Workdata workdata;
    WorkUncdata workunc;
} Args;

/* example: */
/* typedef void (*Intfunc) (gint i); */
/* typedef struct { */
/*     Intfunc f1; */
/* } Fcns */
/* pak delam fcns.f1 = f */


/*
struct {
    gchar *name;

    GtkWidget *gui;
    XXX creategui;
    XXX cleargui;

    ToolData *data;
    XXX initdata;
    XXX cleardata;
    XXX savedata;

    ToolUncData *uncdata;

    XXX saveuncdata;
} Tool;
*/

/*
Tool*
tool_xxx_init(gboolean create_gui) {
initdata;
if (create_gui)
creategui;

}
*/

/*
void
tool_xxx_destroy(Tool *tool)
{
if (tool) {
cleardata;
...
if (tool->gui)
cleargui;
}
}
*/

/*
GSList* tools;
*/

enum {
    R_MODE = 0,
    ER_MODE,
    EIT_MODE,
} Hertzmodes;

#endif

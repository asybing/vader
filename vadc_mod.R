vadc_mod <- '

$PROB Generalized vedotin ADC model

$PLUGIN evtools

$PARAM @annotated
TVCL     : 0.0125/24*70 : Clearance of mAb (L/day/kg)
TVQ      : 0.0247/24*70 : Intercompartmental clearance of mAb (L/day/kg)
TVV1     : 0.0509*70    : Central volume for mAb (L/kg)
TVV2     : 0.0612*70    : Peripheral volume for mAb (L/kg)
TVCLMMAE : 0.861/24*70  : Clearance of unconjugated MMAE (L/day/kg)
TVQMMAE  : 13.7/24*70   : Intercompartmental clearance of uMMAE (L/day/kg)
TVV1MMAE : 0.430*70     : Central volume for uMMAE (L/kg)
TVV2MMAE : 1.24*70      : Peripheral volume for uMMAE (L/kg)
TVKREL   : 0.104/24     : Rate constant for the release of payload from ADC to form naked mAb
TVKDEC   : 0.135/24     : Rate constant of conjugated MMAE releasing payload to form unconjugated MMAE

DAR0   : 4      : Initial drug-antibody ratio
TINF   : 0.5    : Infusion duration (hr)
ADCMW  : 153000 : Molecular weight of BV (g/mol)
MMAEMW : 718    : Molecular weight of vedotin (g/mol)

$CMT @annotated
ADC1       : ADC central compartment. Converted to nanomoles     
ADC2       : ADC peripheral compartment. Converted to nanomoles 
NAB1       : Naked antibody (NAB) central compartment
NAB2       : Naked antibody (NAB) peripheral compartment
ACMMAE1    : acMMAE central compartment. Converted to nanomoles 
ACMMAE2    : acMMAE peripheral compartment. Converted to nanomoles 
MMAE1      : uMMAE central compartment
MMAE2      : uMMAE peripheral compartment
ADC_AUC    : ADC AUC
ACMMAE_AUC : acMMAE AUC
MMAE_AUC   : uMMAE AUC

$EVENT
if(EVID == 1) {

  double ADC_AMT = AMT / ADCMW / 1000 * 1e9; // mg / (g/mol) * (1 g / 1000 mg) * (10^9 nmol / mol)
  
  evt::push(self, evt::infuse(       ADC_AMT, 1,        ADC_AMT/TINF));
  evt::push(self, evt::infuse(DAR0 * ADC_AMT, 5, DAR0 * ADC_AMT/TINF));

}

$MAIN
ADC1_0 = 0;
ADC2_0 = 0;

NAB1_0 = 0;
NAB2_0 = 0;

ACMMAE1_0 = 0;
ACMMAE2_0 = 0;

MMAE1_0 = 0;
MMAE2_0 = 0;

ADC_AUC_0    = 0;
ACMMAE_AUC_0 = 0;
MMAE_AUC_0   = 0;

double CL     = TVCL     * exp(ETA(1));  // Clearance of mAb 
double Q      = TVQ      * exp(ETA(2));  // Intercompartmental clearance of mAb 
double V1     = TVV1     * exp(ETA(3));  // Central compartment of ADC
double V2     = TVV2     * exp(ETA(4));  // Periperal compartment of ADC
double CLMMAE = TVCLMMAE * exp(ETA(5));  // Clearance of uMMAE 
double QMMAE  = TVQMMAE  * exp(ETA(6));  // Intercompartment clearance of uMMAE
double V1MMAE = TVV1MMAE * exp(ETA(7));  // Central compartment of uMMAE
double V2MMAE = TVV2MMAE * exp(ETA(8));  // Peripheral compartment of uMMAE
double KREL   = TVKREL   * exp(ETA(9));  // Rate constant of ADC --> naked mAb
double KDEC   = TVKDEC   * exp(ETA(10)); // Rate constant of acMMAE --> uMMAE

$ODE
double DAR = DAR0;
if((ADC1 + NAB1) > 1e-12) DAR = ACMMAE1/(ADC1 + NAB1);

double CADC1    = ADC1/V1;
double CADC2    = ADC2/V2;
double CNAB1    = NAB1/V1;
double CNAB2    = NAB2/V2;
double CACMMAE1 = ACMMAE1/V1;
double CACMMAE2 = ACMMAE2/V2;
double CMMAE1   = MMAE1/V1MMAE;
double CMMAE2   = MMAE2/V2MMAE;

dxdt_ADC1    = - CL*CADC1 - KREL*ADC1 - Q*CADC1 + Q*CADC2;
dxdt_ADC2    =            - KREL*ADC2 + Q*CADC1 - Q*CADC2;

dxdt_NAB1    = - CL*CNAB1 + KREL*ADC1 - Q*CNAB1 + Q*CNAB2;
dxdt_NAB2    =            + KREL*ADC2 + Q*CNAB1 - Q*CNAB2;

dxdt_ACMMAE1 =                 - DAR*CL*CADC1 - KDEC*ACMMAE1 - Q*CACMMAE1   + Q*CACMMAE2; 
dxdt_ACMMAE2 =                                - KDEC*ACMMAE2 + Q*CACMMAE1   - Q*CACMMAE2; 

dxdt_MMAE1   = - CLMMAE*CMMAE1 + DAR*CL*CADC1 + KDEC*ACMMAE1 - QMMAE*CMMAE1 + QMMAE*CMMAE2;
dxdt_MMAE2   =                                  KDEC*ACMMAE2 + QMMAE*CMMAE1 - QMMAE*CMMAE2;

dxdt_ADC_AUC    = CADC1;
dxdt_ACMMAE_AUC = CACMMAE1;
dxdt_MMAE_AUC   = CMMAE1;

$OMEGA @block @annotated

ECL     : 0.309                                                       : ETA on CL
EQ      : 0     0.307                                                 : ETA on Q
EV1     : 0     0     0.060                                           : ETA on V1
EV2     : 0     0     0     0.534                                     : ETA on V2
ECLMMAE : 0     0     0     0     0.072                               : ETA on CLMMAE
EQMMAE  : 0     0     0     0     0     0                             : ETA on QMMAE
EV1MMAE : 0     0     0     0     0     0     0                       : ETA on V1MMAE
EV2MMAE : 0     0     0     0     0     0     0     0.523             : ETA on V2MMAE
EKREL   : 0     0     0     0     0     0     0     0     0.558       : ETA on KREL
EKDEC   : 0     0     0     0     0     0     0     0     0     0.111 : ETA on KDEC

$CAPTURE 
ADC       =    CADC1*ADCMW*1e-3  // (nmol/L) * (ng / nmol) * (ug / 1000 ng)
ACMMAE    = CACMMAE1*MMAEMW*1e-3
MMAE      =   CMMAE1*MMAEMW*1e-3

AUCADC    =    ADC_AUC*ADCMW*1e-3
AUCACMMAE = ACMMAE_AUC*MMAEMW*1e-3
AUCMMAE   =   MMAE_AUC*MMAEMW*1e-3

DAR0
DAR 
'
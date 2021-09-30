Wed Sep 29 18:12:14 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/TD1/dat47.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m47.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1644.59393080058        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3435E+02 -3.9188E+01 -8.0849E+00 -4.8265E+01  3.7788E+01  6.2019E+01 -3.7322E+00 -1.0799E+00 -2.4023E+01  4.5525E-01
             9.7935E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1651.62807949537        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8423E-01  1.0142E+00  1.0095E+00  1.0594E+00  9.7496E-01  9.2887E-01  1.0144E+00  1.0106E+00  1.0907E+00  9.9464E-01
             9.7586E-01
 PARAMETER:  8.4101E-02  1.1411E-01  1.0949E-01  1.5773E-01  7.4646E-02  2.6218E-02  1.1427E-01  1.1058E-01  1.8683E-01  9.4628E-02
             7.5566E-02
 GRADIENT:   3.0666E+00 -1.7395E+00  4.8490E-01 -8.3460E+00 -3.0110E+00 -1.7723E+00  1.0395E+00 -1.7192E+00 -2.3729E+00  1.5567E+00
            -3.5572E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1652.00156361846        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.8323E-01  8.9441E-01  1.1729E+00  1.1560E+00  9.8677E-01  9.2892E-01  8.7287E-01  1.2375E+00  1.0968E+00  9.9073E-01
             9.6968E-01
 PARAMETER:  8.3089E-02 -1.1588E-02  2.5950E-01  2.4494E-01  8.6686E-02  2.6271E-02 -3.5970E-02  3.1308E-01  1.9236E-01  9.0685E-02
             6.9213E-02
 GRADIENT:   4.5797E+00  1.0795E+01  3.0657E+00  1.3897E+01 -1.2111E+01 -1.1708E+00  1.4441E-01  1.5254E+00  4.0036E+00 -2.9187E+00
            -3.6773E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1652.43118734202        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.8026E-01  8.3533E-01  1.2940E+00  1.1941E+00  1.0246E+00  9.3543E-01  8.3700E-01  1.2517E+00  1.0493E+00  1.0794E+00
             9.7824E-01
 PARAMETER:  8.0058E-02 -7.9927E-02  3.5776E-01  2.7735E-01  1.2431E-01  3.3247E-02 -7.7934E-02  3.2447E-01  1.4808E-01  1.7641E-01
             7.8003E-02
 GRADIENT:  -6.3927E-01  9.5697E+00  3.4693E+00  9.2474E+00 -7.4768E+00  1.7434E+00 -2.4296E-01 -4.4315E-01 -4.0199E+00  9.4299E-01
            -7.6742E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1652.75606272069        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  9.7805E-01  5.9404E-01  1.3925E+00  1.3477E+00  9.8239E-01  9.2842E-01  8.0842E-01  1.2737E+00  9.6874E-01  1.0602E+00
             9.7861E-01
 PARAMETER:  7.7810E-02 -4.2081E-01  4.3109E-01  3.9841E-01  8.2235E-02  2.5730E-02 -1.1268E-01  3.4192E-01  6.8246E-02  1.5847E-01
             7.8382E-02
 GRADIENT:   2.8363E-01  4.1328E+00  5.2123E-01  8.6341E+00 -1.6229E+00 -6.5318E-01 -4.4000E-01  1.2502E-02 -4.2672E-01 -1.7999E-01
            -3.8013E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1652.77754688619        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.7690E-01  5.0870E-01  1.4171E+00  1.4035E+00  9.6310E-01  9.2915E-01  8.4881E-01  1.2762E+00  9.3173E-01  1.0506E+00
             9.7840E-01
 PARAMETER:  7.6625E-02 -5.7590E-01  4.4859E-01  4.3894E-01  6.2399E-02  2.6516E-02 -6.3914E-02  3.4387E-01  2.9287E-02  1.4938E-01
             7.8161E-02
 GRADIENT:  -2.4606E-02  4.3949E+00  1.3457E+00  1.0610E+01 -2.5179E+00 -1.6833E-01 -4.6704E-01 -5.5429E-01 -1.6143E+00 -2.0512E-01
            -5.0438E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1652.79684488696        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  9.7591E-01  4.4048E-01  1.4385E+00  1.4462E+00  9.4965E-01  9.2948E-01  9.1214E-01  1.2911E+00  9.0473E-01  1.0428E+00
             9.7836E-01
 PARAMETER:  7.5620E-02 -7.1990E-01  4.6363E-01  4.6896E-01  4.8334E-02  2.6867E-02  8.0331E-03  3.5551E-01 -1.2062E-04  1.4188E-01
             7.8118E-02
 GRADIENT:  -1.7808E-01  3.5426E+00  1.4209E+00  8.9133E+00 -2.4003E+00  1.1396E-01 -3.9277E-01 -6.8202E-01 -1.8376E+00 -1.6667E-01
            -4.4068E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1652.79900520433        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  9.7542E-01  4.0620E-01  1.4514E+00  1.4678E+00  9.4388E-01  9.2939E-01  9.5679E-01  1.3042E+00  8.9202E-01  1.0390E+00
             9.7833E-01
 PARAMETER:  7.5112E-02 -8.0092E-01  4.7251E-01  4.8376E-01  4.2245E-02  2.6775E-02  5.5826E-02  3.6557E-01 -1.4266E-02  1.3822E-01
             7.8096E-02
 GRADIENT:  -1.8973E-01  3.1077E+00  1.3228E+00  7.9457E+00 -2.2225E+00  1.4836E-01 -3.4113E-01 -6.3860E-01 -1.7181E+00 -1.3588E-01
            -3.8851E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1652.86055789566        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1418
 NPARAMETR:  9.7521E-01  3.9191E-01  1.4527E+00  1.4714E+00  9.4222E-01  9.2892E-01  1.1047E+00  1.3111E+00  8.8780E-01  1.0362E+00
             9.7844E-01
 PARAMETER:  7.4896E-02 -8.3673E-01  4.7341E-01  4.8618E-01  4.0481E-02  2.6272E-02  1.9961E-01  3.7087E-01 -1.9008E-02  1.3556E-01
             7.8201E-02
 GRADIENT:  -7.9355E-02  7.2891E-01 -1.8419E-01 -2.9925E+00  8.4654E-01 -1.7535E-02 -6.4186E-03 -1.0858E-01  2.2232E-01 -3.0118E-02
            -4.2221E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1652.86637637970        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1601
 NPARAMETR:  9.7578E-01  3.9009E-01  1.4520E+00  1.4698E+00  9.4092E-01  9.2917E-01  1.1055E+00  1.3127E+00  8.8672E-01  1.0344E+00
             9.7840E-01
 PARAMETER:  7.5477E-02 -8.4138E-01  4.7296E-01  4.8515E-01  3.9099E-02  2.6540E-02  2.0034E-01  3.7206E-01 -2.0224E-02  1.3383E-01
             7.8168E-02
 GRADIENT:   1.5023E+00  4.8732E-02  7.4365E-02 -7.5506E+00  5.2070E-01  1.0301E-01 -6.3166E-03  2.0327E-02  1.6210E-03 -4.1881E-02
            -2.5836E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1652.86716628153        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1791             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7579E-01  3.9077E-01  1.4503E+00  1.4692E+00  9.4011E-01  9.2918E-01  1.1131E+00  1.3113E+00  8.8701E-01  1.0342E+00
             9.7844E-01
 PARAMETER:  7.5497E-02 -8.3965E-01  4.7175E-01  4.8471E-01  3.8244E-02  2.6544E-02  2.0711E-01  3.7102E-01 -1.9900E-02  1.3361E-01
             7.8204E-02
 GRADIENT:   4.1125E+02  5.0727E+01  9.8027E+00  6.2801E+02  6.9441E+00  3.9534E+01  1.2539E+00  1.9387E+00  1.4975E+01  1.2033E+00
             8.1643E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1652.86759007131        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1981
 NPARAMETR:  9.7580E-01  3.9120E-01  1.4490E+00  1.4688E+00  9.3977E-01  9.2918E-01  1.1148E+00  1.3102E+00  8.8721E-01  1.0338E+00
             9.7843E-01
 PARAMETER:  7.5506E-02 -8.3853E-01  4.7084E-01  4.8447E-01  3.7876E-02  2.6546E-02  2.0863E-01  3.7018E-01 -1.9674E-02  1.3324E-01
             7.8199E-02
 GRADIENT:   1.5259E+00  3.2611E-02  2.2091E-01 -7.8591E+00  6.4021E-02  1.0329E-01  1.8745E-02  5.7589E-02  1.4540E-01  6.7538E-02
             2.1320E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1652.86802430612        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2171             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7581E-01  3.9202E-01  1.4470E+00  1.4684E+00  9.3936E-01  9.2918E-01  1.0996E+00  1.3084E+00  8.8698E-01  1.0329E+00
             9.7835E-01
 PARAMETER:  7.5516E-02 -8.3645E-01  4.6951E-01  4.8419E-01  3.7446E-02  2.6552E-02  1.9493E-01  3.6879E-01 -1.9936E-02  1.3241E-01
             7.8117E-02
 GRADIENT:   4.1127E+02  5.0941E+01  9.6955E+00  6.2663E+02  7.0602E+00  3.9541E+01  1.1222E+00  1.8787E+00  1.4632E+01  1.0920E+00
             7.5741E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1652.86847551193        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2358
 NPARAMETR:  9.7582E-01  3.9232E-01  1.4463E+00  1.4681E+00  9.3911E-01  9.2919E-01  1.1059E+00  1.3077E+00  8.8731E-01  1.0329E+00
             9.7838E-01
 PARAMETER:  7.5524E-02 -8.3568E-01  4.6898E-01  4.8399E-01  3.7173E-02  2.6553E-02  2.0068E-01  3.6830E-01 -1.9557E-02  1.3235E-01
             7.8139E-02
 GRADIENT:   1.5114E+00  9.1269E-02  1.8343E-01 -7.5630E+00  6.5451E-02  1.0387E-01 -9.2461E-03  2.6965E-02 -7.0528E-02  1.1390E-02
            -1.5808E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1652.86873020055        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2543             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7584E-01  3.9264E-01  1.4444E+00  1.4675E+00  9.3886E-01  9.2919E-01  1.1178E+00  1.3067E+00  8.8789E-01  1.0325E+00
             9.7842E-01
 PARAMETER:  7.5543E-02 -8.3487E-01  4.6768E-01  4.8353E-01  3.6911E-02  2.6557E-02  2.1139E-01  3.6753E-01 -1.8903E-02  1.3202E-01
             7.8187E-02
 GRADIENT:   4.1122E+02  5.0718E+01  9.4471E+00  6.2469E+02  7.3790E+00  3.9524E+01  1.3165E+00  1.9511E+00  1.5082E+01  1.1522E+00
             8.4386E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1652.86904207615        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2730             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7584E-01  3.9292E-01  1.4435E+00  1.4674E+00  9.3855E-01  9.2919E-01  1.1171E+00  1.3060E+00  8.8791E-01  1.0323E+00
             9.7842E-01
 PARAMETER:  7.5543E-02 -8.3414E-01  4.6709E-01  4.8350E-01  3.6585E-02  2.6555E-02  2.1074E-01  3.6699E-01 -1.8888E-02  1.3182E-01
             7.8179E-02
 GRADIENT:   4.1119E+02  5.0806E+01  9.4636E+00  6.2476E+02  7.2706E+00  3.9523E+01  1.3084E+00  1.9479E+00  1.5043E+01  1.1645E+00
             8.4093E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1652.86943129543        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2917
 NPARAMETR:  9.7585E-01  3.9342E-01  1.4425E+00  1.4670E+00  9.3827E-01  9.2919E-01  1.1136E+00  1.3048E+00  8.8804E-01  1.0319E+00
             9.7839E-01
 PARAMETER:  7.5557E-02 -8.3287E-01  4.6640E-01  4.8319E-01  3.6286E-02  2.6560E-02  2.0761E-01  3.6608E-01 -1.8744E-02  1.3145E-01
             7.8153E-02
 GRADIENT:   1.5427E+00 -6.6545E-02 -1.6411E-02 -8.1668E+00  3.1360E-01  1.0515E-01  1.6006E-02  6.2345E-02  1.3687E-01  1.6833E-02
             2.4045E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1652.86983934847        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3102             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7586E-01  3.9456E-01  1.4414E+00  1.4666E+00  9.3780E-01  9.2920E-01  1.1062E+00  1.3032E+00  8.8802E-01  1.0316E+00
             9.7833E-01
 PARAMETER:  7.5564E-02 -8.2998E-01  4.6564E-01  4.8295E-01  3.5784E-02  2.6564E-02  2.0097E-01  3.6483E-01 -1.8758E-02  1.3110E-01
             7.8093E-02
 GRADIENT:   4.1122E+02  5.1161E+01  9.7669E+00  6.2351E+02  6.7357E+00  3.9532E+01  1.2039E+00  1.8372E+00  1.4755E+01  1.1282E+00
             7.7522E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1652.87008202577        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3289
 NPARAMETR:  9.7587E-01  3.9476E-01  1.4407E+00  1.4663E+00  9.3775E-01  9.2920E-01  1.1085E+00  1.3029E+00  8.8824E-01  1.0315E+00
             9.7835E-01
 PARAMETER:  7.5572E-02 -8.2949E-01  4.6514E-01  4.8276E-01  3.5732E-02  2.6565E-02  2.0297E-01  3.6456E-01 -1.8516E-02  1.3097E-01
             7.8109E-02
 GRADIENT:   1.5168E+00  7.0245E-02  1.2238E-01 -7.6054E+00  5.5090E-03  1.0368E-01 -1.2385E-03  2.3133E-02 -1.1173E-02  1.3304E-02
            -6.2513E-03

0ITERATION NO.:   91    OBJECTIVE VALUE:  -1652.87008202577        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     3313
 NPARAMETR:  9.7587E-01  3.9476E-01  1.4407E+00  1.4663E+00  9.3775E-01  9.2920E-01  1.1085E+00  1.3029E+00  8.8824E-01  1.0315E+00
             9.7835E-01
 PARAMETER:  7.5572E-02 -8.2949E-01  4.6514E-01  4.8276E-01  3.5732E-02  2.6565E-02  2.0297E-01  3.6456E-01 -1.8516E-02  1.3097E-01
             7.8109E-02
 GRADIENT:  -1.5036E-02 -2.1174E-02  1.5681E-01  3.7576E-01  1.5232E-02 -6.1656E-04 -2.1858E-03  1.4903E-02 -6.2699E-02  1.1258E-02
            -5.9469E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3313
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0268E-04 -1.0069E-02 -3.6377E-02 -2.9209E-03 -4.1645E-02
 SE:             2.9839E-02  6.8312E-03  1.8225E-02  2.8689E-02  2.0756E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9191E-01  1.4049E-01  4.5939E-02  9.1891E-01  4.4808E-02

 ETASHRINKSD(%)  3.5129E-02  7.7115E+01  3.8943E+01  3.8867E+00  3.0465E+01
 ETASHRINKVR(%)  7.0246E-02  9.4763E+01  6.2720E+01  7.6223E+00  5.1649E+01
 EBVSHRINKSD(%)  4.4688E-01  7.7761E+01  4.2763E+01  4.1204E+00  2.7054E+01
 EBVSHRINKVR(%)  8.9177E-01  9.5054E+01  6.7239E+01  8.0710E+00  4.6789E+01
 RELATIVEINF(%)  9.6716E+01  2.1813E-01  8.7032E+00  5.0138E+00  7.8121E+00
 EPSSHRINKSD(%)  4.5749E+01
 EPSSHRINKVR(%)  7.0569E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1652.8700820257707     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -917.71925546203249     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    45.02
 Elapsed covariance  time in seconds:     5.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1652.870       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  3.95E-01  1.44E+00  1.47E+00  9.38E-01  9.29E-01  1.11E+00  1.30E+00  8.88E-01  1.03E+00  9.78E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         2.89E-02  3.99E-01  3.26E-01  2.46E-01  1.64E-01  6.54E-02  1.30E+00  2.93E-01  1.65E-01  2.35E-01  6.47E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        8.36E-04
 
 TH 2
+        4.18E-03  1.59E-01
 
 TH 3
+       -6.80E-04 -2.78E-02  1.06E-01
 
 TH 4
+       -2.48E-03 -9.66E-02  2.18E-02  6.06E-02
 
 TH 5
+        7.65E-04  3.85E-02  3.32E-02 -2.14E-02  2.69E-02
 
 TH 6
+       -7.40E-05 -3.28E-03  3.12E-03  1.91E-03  4.25E-04  4.27E-03
 
 TH 7
+       -1.06E-02 -4.83E-01  6.94E-02  2.89E-01 -1.25E-01  8.72E-03  1.70E+00
 
 TH 8
+       -1.40E-03 -3.67E-02  6.68E-02  2.60E-02  1.36E-02  3.39E-03  1.18E-01  8.58E-02
 
 TH 9
+        1.82E-03  5.83E-02 -1.60E-02 -3.53E-02  1.21E-02 -1.19E-03 -1.67E-01 -1.73E-02  2.73E-02
 
 TH10
+        5.92E-04  4.42E-02  4.31E-02 -2.48E-02  3.32E-02  2.24E-03 -1.44E-01  1.13E-02  1.31E-02  5.54E-02
 
 TH11
+        1.09E-04  6.78E-04  3.25E-03 -8.05E-05  1.47E-03  1.89E-04 -4.57E-03  3.17E-04 -2.97E-04  1.50E-03  4.19E-03
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.89E-02
 
 TH 2
+        3.62E-01  3.99E-01
 
 TH 3
+       -7.22E-02 -2.14E-01  3.26E-01
 
 TH 4
+       -3.49E-01 -9.85E-01  2.71E-01  2.46E-01
 
 TH 5
+        1.61E-01  5.88E-01  6.20E-01 -5.30E-01  1.64E-01
 
 TH 6
+       -3.91E-02 -1.26E-01  1.46E-01  1.19E-01  3.96E-02  6.54E-02
 
 TH 7
+       -2.81E-01 -9.29E-01  1.63E-01  9.00E-01 -5.84E-01  1.02E-01  1.30E+00
 
 TH 8
+       -1.65E-01 -3.15E-01  6.99E-01  3.60E-01  2.84E-01  1.77E-01  3.09E-01  2.93E-01
 
 TH 9
+        3.81E-01  8.86E-01 -2.97E-01 -8.69E-01  4.48E-01 -1.10E-01 -7.74E-01 -3.57E-01  1.65E-01
 
 TH10
+        8.71E-02  4.71E-01  5.62E-01 -4.28E-01  8.61E-01  1.46E-01 -4.69E-01  1.65E-01  3.38E-01  2.35E-01
 
 TH11
+        5.84E-02  2.63E-02  1.54E-01 -5.06E-03  1.39E-01  4.47E-02 -5.41E-02  1.67E-02 -2.78E-02  9.84E-02  6.47E-02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.59E+03
 
 TH 2
+       -2.05E+02  5.47E+02
 
 TH 3
+       -1.21E+02  6.79E+01  1.24E+02
 
 TH 4
+       -1.08E+02  5.47E+02 -2.65E+00  7.39E+02
 
 TH 5
+        1.90E+02 -1.81E+02 -2.27E+02 -4.92E+01  7.00E+02
 
 TH 6
+       -5.25E+01  6.92E+01  1.96E+01  5.52E+01  7.00E+00  2.70E+02
 
 TH 7
+       -1.50E+01  3.13E+01  1.02E-01  1.87E+01  5.25E+00  4.16E+00  5.52E+00
 
 TH 8
+        4.05E+01 -3.27E+01 -2.61E+01 -2.64E+01 -1.69E-01 -1.95E+01 -3.08E+00  3.11E+01
 
 TH 9
+       -5.46E+01 -1.47E+02  2.48E+01 -8.69E+01 -5.55E+01 -2.05E+01 -1.10E+01  9.26E+00  2.15E+02
 
 TH10
+        4.66E+01 -1.22E+01 -1.49E+01  3.09E-01 -9.60E+01 -4.01E+01 -2.26E+00  1.81E+01  1.22E+01  8.54E+01
 
 TH11
+       -2.09E+01 -3.44E+01 -1.57E+01 -3.86E+01 -1.11E+01 -1.98E+01 -1.60E-01  1.34E+01  2.31E+01  1.45E+01  2.57E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       50.876
Stop Time:
Wed Sep 29 18:13:06 CDT 2021

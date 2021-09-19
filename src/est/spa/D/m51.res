Sat Sep 18 15:25:27 CDT 2021
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
$DATA ../../../../data/spa/D/dat51.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m51.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12451.6463446620        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.1769E+02  2.5015E+02  1.1808E+01  1.2098E+02  5.1309E+01 -1.6004E+03 -8.1338E+02 -7.7364E+01 -1.1433E+03 -4.0954E+02
            -2.4076E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -605.386817327995        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2528E+00  1.0808E+00  9.5308E-01  1.5068E+00  1.2424E+00  1.9367E+00  1.2363E+00  9.9603E-01  1.1502E+00  1.0376E+00
             1.4511E+01
 PARAMETER:  3.2542E-01  1.7770E-01  5.1945E-02  5.1002E-01  3.1703E-01  7.6099E-01  3.1213E-01  9.6022E-02  2.3993E-01  1.3689E-01
             2.7749E+00
 GRADIENT:  -2.3231E+01  6.3869E+00 -7.0749E+00  6.9733E+00 -3.1034E+00  5.1216E+01  2.9713E+00  4.3013E+00  1.2737E+01  3.0028E+00
             1.6292E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -614.036332634453        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.2818E+00  8.4881E-01  1.0769E+00  1.8416E+00  2.2192E+00  1.7115E+00  3.9966E+00  7.1548E-01  9.6372E-01  1.9216E+00
             1.3588E+01
 PARAMETER:  3.4825E-01 -6.3915E-02  1.7405E-01  7.1065E-01  8.9715E-01  6.3736E-01  1.4854E+00 -2.3480E-01  6.3043E-02  7.5317E-01
             2.7092E+00
 GRADIENT:   1.2320E+01  3.9721E+01 -8.9903E+00  2.7195E+01 -1.0355E+01  1.6867E+01  3.3216E+01  1.7507E+00  3.3862E+00  2.5643E-01
             8.9465E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -642.847892765639        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.1097E+00  3.2458E-01  1.4502E+00  1.6495E+00  1.8007E+00  1.5495E+00  4.4277E+00  7.4978E-01  6.0501E-01  2.3604E+00
             1.1528E+01
 PARAMETER:  2.0405E-01 -1.0252E+00  4.7169E-01  6.0047E-01  6.8816E-01  5.3790E-01  1.5879E+00 -1.8797E-01 -4.0251E-01  9.5885E-01
             2.5448E+00
 GRADIENT:  -9.7020E+00  1.1215E+01  1.5655E+01 -3.7929E+01 -1.2796E+01  1.0148E+01  1.4413E+01  6.5831E-01  6.4098E-01  4.4379E-01
             3.8270E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -674.166628401115        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.5770E-01  1.2623E-01  7.3409E-01  1.4581E+00  1.3572E+01  1.6173E+00  6.3468E+00  1.2576E-02  1.2054E-01  1.1442E+01
             1.0950E+01
 PARAMETER:  5.6774E-02 -1.9696E+00 -2.0912E-01  4.7715E-01  2.7080E+00  5.8075E-01  1.9480E+00 -4.2760E+00 -2.0158E+00  2.5373E+00
             2.4933E+00
 GRADIENT:  -2.5738E+01  6.1594E+00  5.3185E+01 -9.2928E+01  4.2046E+00  2.7096E+01  1.4751E+01  2.3987E-04  5.2833E-01 -2.3994E+00
             2.4495E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -684.384448122697        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      381
 NPARAMETR:  9.2687E-01  8.6354E-02  4.9035E-01  1.4016E+00  2.3850E+01  1.6402E+00  6.9239E+00  1.0000E-02  6.4363E-02  1.4617E+01
             1.0996E+01
 PARAMETER:  2.4054E-02 -2.3493E+00 -6.1264E-01  4.3760E-01  3.2718E+00  5.9481E-01  2.0350E+00 -6.3924E+00 -2.6432E+00  2.7822E+00
             2.4975E+00
 GRADIENT:  -2.2805E+01  3.9294E+00  1.5962E+01 -9.1080E+00  6.8403E-01  1.5124E+01  4.0047E+00  0.0000E+00  1.2328E-01  1.4671E+01
             1.6965E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -684.814603499458        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      464
 NPARAMETR:  9.2786E-01  8.5102E-02  4.8591E-01  1.3996E+00  2.3636E+01  1.6374E+00  6.9238E+00  1.0000E-02  6.4059E-02  1.4388E+01
             1.0978E+01
 PARAMETER:  2.5129E-02 -2.3639E+00 -6.2173E-01  4.3620E-01  3.2628E+00  5.9313E-01  2.0350E+00 -6.4292E+00 -2.6479E+00  2.7664E+00
             2.4959E+00
 GRADIENT:  -2.1227E+01  3.8565E+00  1.5036E+01 -7.2924E+00  7.3352E-01  1.4469E+01  3.8351E+00  0.0000E+00  1.1855E-01  1.4499E+01
             1.4972E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -685.235679170852        NO. OF FUNC. EVALS.: 153
 CUMULATIVE NO. OF FUNC. EVALS.:      617
 NPARAMETR:  9.2808E-01  8.5072E-02  4.8442E-01  1.4010E+00  2.3446E+01  1.6328E+00  6.8599E+00  1.0000E-02  1.0000E-02  1.4154E+01
             1.0654E+01
 PARAMETER:  2.5367E-02 -2.3643E+00 -6.2480E-01  4.3717E-01  3.2547E+00  5.9031E-01  2.0257E+00 -6.4292E+00 -9.4809E+00  2.7500E+00
             2.4660E+00
 GRADIENT:  -2.0483E+01  3.9639E+00  1.2632E+01 -7.4610E+00  1.6056E+00  8.7895E+00  3.7215E+00  0.0000E+00  0.0000E+00  5.9478E+00
            -8.0064E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -688.192409348935        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      784
 NPARAMETR:  9.5096E-01  1.0000E-02  4.5342E-01  1.4013E+00  2.3391E+01  1.5844E+00  3.9287E+00  1.0000E-02  1.0000E-02  1.4176E+01
             1.0706E+01
 PARAMETER:  4.9714E-02 -4.7618E+00 -6.9093E-01  4.3741E-01  3.2524E+00  5.6019E-01  1.4683E+00 -6.4292E+00 -8.9392E+00  2.7516E+00
             2.4708E+00
 GRADIENT:   9.1764E+00  0.0000E+00  1.7493E-01 -1.3693E+01  1.7698E+00  1.2278E+00  6.6813E-03  0.0000E+00  0.0000E+00  1.1215E+01
            -2.7588E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -688.380189945839        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:      915             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3760E-01  1.0000E-02  4.5446E-01  1.4004E+00  2.3340E+01  1.5633E+00  3.8368E+00  1.0000E-02  1.0000E-02  1.4056E+01
             1.0762E+01
 PARAMETER:  3.5566E-02 -4.8119E+00 -6.8864E-01  4.3675E-01  3.2502E+00  5.4677E-01  1.4446E+00 -6.4292E+00 -8.9392E+00  2.7431E+00
             2.4761E+00
 GRADIENT:   1.3093E+00  0.0000E+00  2.5412E-01  9.4512E+00  5.8340E-01 -3.5711E-01  6.1393E-03  0.0000E+00  0.0000E+00  1.5858E+01
            -8.2554E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -688.401820565089        NO. OF FUNC. EVALS.: 113
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  9.3590E-01  1.0000E-02  4.5226E-01  1.3958E+00  2.3331E+01  1.5634E+00  3.8207E+00  1.0000E-02  1.0000E-02  1.4059E+01
             1.0760E+01
 PARAMETER:  3.3746E-02 -4.8119E+00 -6.9281E-01  4.3340E-01  3.2501E+00  5.4691E-01  1.4413E+00 -6.4292E+00 -8.9392E+00  2.7431E+00
             2.4760E+00
 GRADIENT:  -4.4432E+03  0.0000E+00  5.2464E-01 -1.0155E+03  6.8702E+01  8.1203E+02  5.5930E-03  0.0000E+00  0.0000E+00 -1.5055E+02
             1.6876E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1028
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4812E-02 -8.4231E-04  2.0655E-05 -6.0762E-04 -1.6516E-02
 SE:             2.7388E-02  2.9667E-04  1.1112E-04  3.3495E-04  5.9600E-03
 N:                     100         100         100         100         100

 P VAL.:         5.8862E-01  4.5223E-03  8.5253E-01  6.9668E-02  5.5875E-03

 ETASHRINKSD(%)  8.2462E+00  9.9006E+01  9.9628E+01  9.8878E+01  8.0033E+01
 ETASHRINKVR(%)  1.5812E+01  9.9990E+01  9.9999E+01  9.9987E+01  9.6013E+01
 EBVSHRINKSD(%)  4.9090E+00  9.8960E+01  9.9638E+01  9.8926E+01  7.8720E+01
 EBVSHRINKVR(%)  9.5770E+00  9.9989E+01  9.9999E+01  9.9988E+01  9.5472E+01
 RELATIVEINF(%)  4.3859E+01  2.5087E-04  4.9352E-05  1.1140E-04  3.7351E+00
 EPSSHRINKSD(%)  4.6282E+00
 EPSSHRINKVR(%)  9.0422E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -688.40182056508877     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       46.749005998649409     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.78
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -688.402       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.36E-01  1.00E-02  4.53E-01  1.40E+00  2.33E+01  1.56E+00  3.82E+00  1.00E-02  1.00E-02  1.41E+01  1.08E+01
 


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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.27E+07
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -3.79E+06  0.00E+00  1.14E+06
 
 TH 4
+        9.54E+02  0.00E+00 -3.05E+03  3.00E+05
 
 TH 5
+       -8.94E+00  0.00E+00  2.16E+01  2.89E+01  4.52E+01
 
 TH 6
+       -7.27E+02  0.00E+00  1.96E+03 -3.40E+02  2.55E+00  1.52E+05
 
 TH 7
+        1.08E+00  0.00E+00 -7.19E-02 -2.04E-01 -1.36E-03 -1.15E-01  1.47E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        1.77E+01  0.00E+00 -4.44E+01 -5.79E+01 -5.14E+01 -5.10E+00  3.86E-03  0.00E+00  0.00E+00  1.79E+02
 
 TH11
+       -3.44E+01  0.00E+00  7.00E+01  6.52E+01  7.17E+01  9.08E+00  2.22E-03  0.00E+00  0.00E+00 -1.42E+02  3.63E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       22.963
Stop Time:
Sat Sep 18 15:25:52 CDT 2021

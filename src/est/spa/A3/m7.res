Sat Sep 25 09:02:51 CDT 2021
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
$DATA ../../../../data/spa/A3/dat7.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   347.475135436739        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2558E+02  9.8282E+01  1.0364E+02  9.7953E+00  1.9579E+02 -4.1511E+01 -6.5538E+01 -7.2531E+01 -1.3971E+02 -1.9844E+02
            -3.4541E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -946.144582977227        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1414E+00  9.1711E-01  7.0737E-01  1.4371E+00  6.5655E-01  9.0139E-01  1.0738E+00  1.2222E+00  1.2618E+00  1.4068E+00
             1.0203E+01
 PARAMETER:  2.3229E-01  1.3472E-02 -2.4620E-01  4.6261E-01 -3.2075E-01 -3.8154E-03  1.7121E-01  3.0068E-01  3.3252E-01  4.4134E-01
             2.4226E+00
 GRADIENT:  -6.6035E+01  8.0994E+00 -1.6865E+01  5.0167E+01 -2.2199E+01 -1.4734E+01  8.5598E+00  1.0363E+01  3.8788E+01  3.4996E+01
             3.7198E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1018.94042900884        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.6264E-01  5.5207E-01  1.3059E-01  1.2390E+00  1.9752E-01  1.4603E+00  1.3736E+00  1.6437E+00  2.2316E+00  7.2982E-01
             6.5788E+00
 PARAMETER:  6.1925E-02 -4.9408E-01 -1.9357E+00  3.1432E-01 -1.5219E+00  4.7865E-01  4.1741E-01  5.9698E-01  9.0271E-01 -2.1496E-01
             1.9839E+00
 GRADIENT:  -9.9632E+01  1.0548E+02  5.3835E+01  7.0545E+01 -1.5934E+02  5.9327E+01  2.2595E+01  1.2207E+01  4.9711E+01  1.9750E+01
             2.4882E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1200.70181645761        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  8.8687E-01  5.5044E-01  9.6583E-02  7.6346E-01  2.3228E-01  1.2305E+00  2.3812E-01  1.0858E+00  1.7555E+00  3.8677E-01
             3.4828E+00
 PARAMETER: -2.0062E-02 -4.9704E-01 -2.2374E+00 -1.6990E-01 -1.3598E+00  3.0744E-01 -1.3350E+00  1.8231E-01  6.6274E-01 -8.4992E-01
             1.3478E+00
 GRADIENT:  -1.0768E+02  3.4185E+01 -2.1705E+01 -6.8463E+01 -2.9328E+01  2.5335E+01  6.5164E-01 -6.3522E+00  6.4782E+00 -9.5846E+00
             3.0929E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1226.50393976979        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.7244E-01  5.1759E-01  1.3081E-01  9.2219E-01  2.4553E-01  1.0899E+00  2.4362E-01  1.4599E+00  1.4925E+00  4.8051E-01
             3.0697E+00
 PARAMETER:  7.2049E-02 -5.5857E-01 -1.9340E+00  1.8991E-02 -1.3043E+00  1.8607E-01 -1.3121E+00  4.7834E-01  5.0042E-01 -6.3290E-01
             1.2216E+00
 GRADIENT:   1.4324E+00 -4.0713E+00 -1.7942E+00 -3.8081E+00  1.1299E+01 -6.8603E-02  9.2982E-01 -2.5029E+00 -2.8298E+00  7.8167E-02
             2.8307E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1226.60144596151        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.6935E-01  5.0508E-01  1.2728E-01  9.2657E-01  2.3900E-01  1.0922E+00  1.6227E-01  1.5055E+00  1.5356E+00  4.9244E-01
             3.0393E+00
 PARAMETER:  6.8867E-02 -5.8305E-01 -1.9613E+00  2.3729E-02 -1.3313E+00  1.8823E-01 -1.7185E+00  5.0914E-01  5.2892E-01 -6.0838E-01
             1.2116E+00
 GRADIENT:  -2.9008E+00  3.6431E+00  1.9841E-01  3.0096E+00 -3.1567E+00  6.8959E-01  4.3026E-01 -1.3768E-02  3.4416E-01 -8.5082E-01
            -8.3260E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1227.17647116972        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.7321E-01  5.4973E-01  1.3771E-01  9.2886E-01  2.5866E-01  1.0894E+00  8.9921E-02  1.5448E+00  1.5126E+00  4.6887E-01
             3.0770E+00
 PARAMETER:  7.2849E-02 -4.9833E-01 -1.8826E+00  2.6198E-02 -1.2522E+00  1.8560E-01 -2.3088E+00  5.3490E-01  5.1386E-01 -6.5744E-01
             1.2239E+00
 GRADIENT:   1.1078E+00 -1.6375E+00 -4.4074E-01 -1.6295E+00  3.1807E+00  6.3076E-01  1.0744E-01  5.8524E-02  1.2487E-01  5.4472E-01
             1.1547E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1227.23516290925        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.7278E-01  5.4640E-01  1.3734E-01  9.3101E-01  2.5731E-01  1.0882E+00  1.1672E-02  1.5456E+00  1.5153E+00  4.6825E-01
             3.0691E+00
 PARAMETER:  7.2404E-02 -5.0441E-01 -1.8853E+00  2.8516E-02 -1.2575E+00  1.8452E-01 -4.3505E+00  5.3543E-01  5.1564E-01 -6.5876E-01
             1.2214E+00
 GRADIENT:   2.2871E-01  3.3595E-01  1.4706E-02  1.8954E-01 -4.0640E-01  1.3778E-01  1.7433E-03  5.3279E-02 -7.5979E-02  9.4316E-03
            -7.5438E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1227.23551359410        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  9.7270E-01  5.4612E-01  1.3736E-01  9.3095E-01  2.5726E-01  1.0879E+00  1.0000E-02  1.5443E+00  1.5153E+00  4.6816E-01
             3.0694E+00
 PARAMETER:  7.2316E-02 -5.0492E-01 -1.8851E+00  2.8449E-02 -1.2577E+00  1.8424E-01 -4.6206E+00  5.3457E-01  5.1560E-01 -6.5894E-01
             1.2215E+00
 GRADIENT:   2.8932E-02  7.0164E-02  2.5991E-02  5.3268E-02 -1.1978E-01  3.1170E-02  0.0000E+00 -1.1621E-02 -6.1730E-02  8.6421E-03
            -3.4845E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      762
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2829E-03 -4.0105E-04  1.8681E-02 -4.8062E-03  1.8652E-02
 SE:             2.8826E-02  1.9940E-04  1.9250E-02  2.6131E-02  1.5949E-02
 N:                     100         100         100         100         100

 P VAL.:         9.0933E-01  4.4299E-02  3.3184E-01  8.5407E-01  2.4220E-01

 ETASHRINKSD(%)  3.4279E+00  9.9332E+01  3.5509E+01  1.2459E+01  4.6569E+01
 ETASHRINKVR(%)  6.7383E+00  9.9996E+01  5.8409E+01  2.3365E+01  7.1452E+01
 EBVSHRINKSD(%)  3.1688E+00  9.9282E+01  3.5600E+01  1.1069E+01  4.6859E+01
 EBVSHRINKVR(%)  6.2372E+00  9.9995E+01  5.8527E+01  2.0912E+01  7.1760E+01
 RELATIVEINF(%)  9.0264E+01  2.4619E-04  1.2088E+01  5.8240E+01  1.0313E+00
 EPSSHRINKSD(%)  3.4121E+01
 EPSSHRINKVR(%)  5.6600E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1227.2355135940952     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -492.08468703035703     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1227.236       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.73E-01  5.46E-01  1.37E-01  9.31E-01  2.57E-01  1.09E+00  1.00E-02  1.54E+00  1.52E+00  4.68E-01  3.07E+00
 


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
+        9.36E+02
 
 TH 2
+        1.63E+00  2.15E+03
 
 TH 3
+       -2.96E+02  2.63E+03  1.03E+04
 
 TH 4
+       -2.47E+01  2.08E+02 -4.20E+02  4.43E+02
 
 TH 5
+        1.90E+02 -6.55E+03 -1.07E+04 -3.70E+02  2.25E+04
 
 TH 6
+        2.67E+00 -1.24E+01  9.53E+00 -9.89E+00  5.74E+01  1.47E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        2.24E+00 -4.73E+00 -1.51E+01  1.78E-01 -1.87E+01  3.72E+00  0.00E+00  1.72E+01
 
 TH 9
+        9.77E+00 -4.49E+01  7.48E+01 -7.80E+00  1.15E+02  2.01E+00  0.00E+00  1.01E+00  4.54E+01
 
 TH10
+       -5.06E-01 -1.10E+02  2.50E+01  9.19E-01  4.18E+02  5.01E+00  0.00E+00  3.87E+00  9.08E+00  1.10E+02
 
 TH11
+       -1.65E+01 -3.16E+01 -5.77E+00 -2.37E+00  5.69E+01 -1.80E-01  0.00E+00  5.93E+00  7.22E+00  1.43E+01  2.79E+01
 
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
 #CPUT: Total CPU Time in Seconds,       15.704
Stop Time:
Sat Sep 25 09:03:08 CDT 2021

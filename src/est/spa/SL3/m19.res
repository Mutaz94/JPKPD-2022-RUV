Wed Sep 29 16:25:23 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat19.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1688.67863569427        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9391E+02 -2.0868E+01 -3.2759E+01  6.5487E+01  1.0457E+02  5.7989E+01  6.0489E+00  1.0760E-01  4.5949E+01 -2.7545E+01
            -1.9932E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1698.44141991868        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0478E+00  1.0747E+00  9.4672E-01  9.7551E-01  9.6608E-01  9.6235E-01  9.6566E-01  1.0078E+00  7.4700E-01  1.1293E+00
             1.0786E+00
 PARAMETER:  1.4673E-01  1.7208E-01  4.5249E-02  7.5207E-02  6.5496E-02  6.1626E-02  6.5053E-02  1.0776E-01 -1.9169E-01  2.2158E-01
             1.7563E-01
 GRADIENT:   3.8274E+01 -1.0145E+01 -1.4295E+01  7.5601E+00  2.0209E+01 -5.8452E+00 -8.4667E+00  3.3166E+00 -8.3516E+00  4.3575E+00
             1.1359E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1699.06778573314        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0414E+00  1.0179E+00  8.6961E-01  1.0108E+00  9.0130E-01  9.7786E-01  1.1023E+00  8.6228E-01  6.9918E-01  1.0603E+00
             1.0587E+00
 PARAMETER:  1.4059E-01  1.1773E-01 -3.9708E-02  1.1078E-01 -3.9200E-03  7.7615E-02  1.9741E-01 -4.8177E-02 -2.5785E-01  1.5851E-01
             1.5706E-01
 GRADIENT:   2.1822E+01  2.1811E+00 -1.9864E+01  2.4887E+01  2.4359E+01  6.0198E-01 -2.1379E+00  3.5473E+00 -8.6710E+00  3.0152E+00
             5.0001E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.48934819894        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.0320E+00  9.5234E-01  7.6117E-01  1.0317E+00  8.0029E-01  9.7637E-01  1.1532E+00  4.8090E-01  7.4368E-01  9.8331E-01
             1.0412E+00
 PARAMETER:  1.3147E-01  5.1170E-02 -1.7290E-01  1.3116E-01 -1.2278E-01  7.6081E-02  2.4250E-01 -6.3209E-01 -1.9614E-01  8.3172E-02
             1.4041E-01
 GRADIENT:  -1.1117E+00 -3.9419E-01 -2.4827E+00  5.1183E-01  1.0591E-01 -1.2227E-01  5.1137E-01  8.3086E-01  1.1781E+00  1.2208E+00
             1.1439E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1700.58340681038        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      728
 NPARAMETR:  1.0320E+00  8.6451E-01  7.4573E-01  1.0783E+00  7.5713E-01  9.7589E-01  1.2510E+00  3.1005E-01  7.1438E-01  9.6456E-01
             1.0410E+00
 PARAMETER:  1.3148E-01 -4.5596E-02 -1.9339E-01  1.7540E-01 -1.7822E-01  7.5599E-02  3.2398E-01 -1.0710E+00 -2.3634E-01  6.3913E-02
             1.4022E-01
 GRADIENT:  -8.4186E-01 -1.7989E+00  3.8968E-01 -3.0511E+00 -5.7116E-01 -5.8922E-02 -7.5527E-02  1.0443E-01  1.7058E-02  2.0339E-01
             1.9211E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1700.72044819574        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      906
 NPARAMETR:  1.0328E+00  1.0011E+00  6.8408E-01  9.9551E-01  7.7865E-01  9.7752E-01  1.0986E+00  9.4013E-02  7.6138E-01  9.6543E-01
             1.0408E+00
 PARAMETER:  1.3229E-01  1.0113E-01 -2.7969E-01  9.5502E-02 -1.5020E-01  7.7266E-02  1.9402E-01 -2.2643E+00 -1.7262E-01  6.4820E-02
             1.3994E-01
 GRADIENT:  -1.5504E+00  6.4096E-01 -8.0379E-01  1.8549E+00  1.1219E-01 -1.3078E-01 -1.6959E-01  2.3492E-02 -2.8071E-02  4.7846E-01
             1.2644E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1700.73345807607        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1080
 NPARAMETR:  1.0348E+00  1.0194E+00  6.8006E-01  9.8367E-01  7.8427E-01  9.7825E-01  1.0833E+00  4.8056E-02  7.6818E-01  9.6506E-01
             1.0407E+00
 PARAMETER:  1.3419E-01  1.1921E-01 -2.8558E-01  8.3532E-02 -1.4300E-01  7.8006E-02  1.7997E-01 -2.9354E+00 -1.6373E-01  6.4433E-02
             1.3991E-01
 GRADIENT:   2.8490E+00  2.3085E-01  6.2216E-01 -6.3870E-01 -6.9211E-01  1.2516E-01 -3.2542E-02  5.2155E-03 -6.2364E-02 -1.8301E-01
            -1.1541E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1700.73602429208        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1255
 NPARAMETR:  1.0347E+00  1.0188E+00  6.8026E-01  9.8436E-01  7.8449E-01  9.7841E-01  1.0840E+00  1.2473E-02  7.6842E-01  9.6674E-01
             1.0409E+00
 PARAMETER:  1.3409E-01  1.1859E-01 -2.8529E-01  8.4238E-02 -1.4272E-01  7.8176E-02  1.8064E-01 -4.2842E+00 -1.6342E-01  6.6178E-02
             1.4004E-01
 GRADIENT:   2.5898E+00  2.5333E-01  3.1688E-02  2.9601E-01  7.5900E-02  1.9240E-01  2.9642E-02  3.9524E-04  5.7826E-02  2.5836E-02
             5.2895E-03

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1700.73606335009        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1354
 NPARAMETR:  1.0353E+00  1.0171E+00  6.8030E-01  9.8374E-01  7.8406E-01  9.7750E-01  1.0837E+00  1.0000E-02  7.6767E-01  9.6654E-01
             1.0407E+00
 PARAMETER:  1.3405E-01  1.1809E-01 -2.8493E-01  8.4539E-02 -1.4277E-01  7.8148E-02  1.8115E-01 -4.8469E+00 -1.6372E-01  6.6185E-02
             1.4007E-01
 GRADIENT:  -3.0899E-01  2.3802E-01  4.0611E-02  4.2073E-01  1.6055E-01  7.0896E-02  1.8500E-02  0.0000E+00  2.2104E-02  8.0981E-03
             1.3698E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1354
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.7424E-06 -4.6453E-03 -4.1545E-04 -1.6874E-03 -1.3279E-02
 SE:             2.9808E-02  2.1407E-02  1.7980E-04  2.3327E-02  2.4202E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9977E-01  8.2821E-01  2.0852E-02  9.4233E-01  5.8322E-01

 ETASHRINKSD(%)  1.4011E-01  2.8284E+01  9.9398E+01  2.1852E+01  1.8920E+01
 ETASHRINKVR(%)  2.8003E-01  4.8568E+01  9.9996E+01  3.8929E+01  3.4261E+01
 EBVSHRINKSD(%)  4.8270E-01  2.8237E+01  9.9474E+01  2.2374E+01  1.7422E+01
 EBVSHRINKVR(%)  9.6306E-01  4.8501E+01  9.9997E+01  3.9742E+01  3.1808E+01
 RELATIVEINF(%)  9.8712E+01  2.5057E+00  2.9279E-04  3.1366E+00  6.2521E+00
 EPSSHRINKSD(%)  4.3812E+01
 EPSSHRINKVR(%)  6.8429E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1700.7360633500869     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -965.58523678634867     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.09
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1700.736       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.02E+00  6.80E-01  9.85E-01  7.84E-01  9.78E-01  1.08E+00  1.00E-02  7.68E-01  9.67E-01  1.04E+00
 


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
+        1.07E+03
 
 TH 2
+       -8.70E+00  4.88E+02
 
 TH 3
+        1.93E+01  2.04E+02  7.55E+02
 
 TH 4
+       -1.84E+01  4.76E+02 -3.91E+02  1.17E+03
 
 TH 5
+       -3.11E+00 -3.72E+02 -8.44E+02  3.83E+02  1.31E+03
 
 TH 6
+        3.13E-01 -2.06E+00  5.55E+00 -4.51E+00 -1.50E+00  2.04E+02
 
 TH 7
+        1.42E+00  3.01E+01 -1.32E+01 -9.53E+00 -3.05E+00  2.93E-01  5.19E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.99E+00 -2.37E+01 -2.79E+01  2.36E+01  9.01E+00 -5.40E-01  2.51E+01  0.00E+00  1.40E+02
 
 TH10
+       -1.05E+00 -1.18E+01 -7.20E+01 -1.36E+01 -4.96E+01 -1.02E-01  1.26E+01  0.00E+00  9.57E+00  9.88E+01
 
 TH11
+       -7.89E+00 -1.55E+01 -3.58E+01 -2.39E+00  9.53E+00  2.23E+00  5.58E+00  0.00E+00  1.32E+01  1.82E+01  1.98E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.725
Stop Time:
Wed Sep 29 16:25:48 CDT 2021

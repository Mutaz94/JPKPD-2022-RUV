Sat Sep 18 11:43:49 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat50.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1692.32455021990        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2751E+01 -3.8997E+01 -4.1579E+01 -3.9584E+01  4.4611E+00  4.8043E+00  1.5159E+01  1.9134E+01 -1.8572E+00  4.2820E+01
            -8.0582E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1704.16431467687        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0142E+00  9.5687E-01  1.1888E+00  1.0218E+00  1.0492E+00  9.8336E-01  8.1846E-01  8.5042E-01  1.0504E+00  7.5170E-01
             1.0724E+00
 PARAMETER:  1.1407E-01  5.5915E-02  2.7292E-01  1.2161E-01  1.4802E-01  8.3223E-02 -1.0034E-01 -6.2020E-02  1.4919E-01 -1.8542E-01
             1.6990E-01
 GRADIENT:   8.7771E+01 -3.8616E+01 -9.2362E-01 -5.5792E+01  2.7308E+01 -2.6352E+00  7.9040E+00  4.7173E+00  1.9827E+00 -4.2851E+00
             1.0632E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1707.12510485423        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0093E+00  8.1310E-01  1.1475E+00  1.1476E+00  9.9280E-01  9.7909E-01  6.9268E-01  3.7940E-01  9.6349E-01  8.8540E-01
             1.0054E+00
 PARAMETER:  1.0925E-01 -1.0691E-01  2.3760E-01  2.3770E-01  9.2778E-02  7.8870E-02 -2.6719E-01 -8.6916E-01  6.2805E-02 -2.1719E-02
             1.0537E-01
 GRADIENT:   8.2502E+01 -8.1619E+00 -1.9021E+01  1.5558E+01  3.8009E+01 -3.3756E+00  1.2035E+00  8.5266E-01 -6.9253E+00  9.0647E+00
            -1.9013E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1708.91307636016        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.8839E-01  8.1930E-01  1.0550E+00  1.1390E+00  9.3277E-01  9.9467E-01  7.4685E-01  3.5469E-01  9.7846E-01  7.8198E-01
             1.0370E+00
 PARAMETER:  8.8319E-02 -9.9301E-02  1.5353E-01  2.3017E-01  3.0406E-02  9.4657E-02 -1.9189E-01 -9.3651E-01  7.8229E-02 -1.4593E-01
             1.3631E-01
 GRADIENT:   2.6445E+01 -3.9123E+00 -1.0743E+01  1.3676E+01  1.6318E+01  3.1418E+00  1.0093E+00  6.3945E-01 -4.7070E-01  1.4521E+00
            -4.6076E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1709.36791628505        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      407
 NPARAMETR:  9.9199E-01  7.0774E-01  1.1001E+00  1.2112E+00  9.0809E-01  9.9534E-01  6.3116E-01  1.7797E-01  9.4773E-01  8.0361E-01
             1.0551E+00
 PARAMETER:  9.1958E-02 -2.4567E-01  1.9544E-01  2.9158E-01  3.5890E-03  9.5333E-02 -3.6019E-01 -1.6261E+00  4.6309E-02 -1.1864E-01
             1.5359E-01
 GRADIENT:  -1.2867E+00 -2.9649E-01 -1.5914E+00  4.3648E-01  1.4210E+00  2.1460E-01 -5.2401E-01  1.7223E-01  4.1291E-01 -1.4835E-01
            -1.4817E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1709.50187026360        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      583
 NPARAMETR:  9.9102E-01  5.8699E-01  1.0938E+00  1.2851E+00  8.6258E-01  9.9362E-01  7.8885E-01  5.0782E-02  8.8912E-01  7.9499E-01
             1.0531E+00
 PARAMETER:  9.0976E-02 -4.3275E-01  1.8970E-01  3.5086E-01 -4.7822E-02  9.3596E-02 -1.3718E-01 -2.8802E+00 -1.7518E-02 -1.2942E-01
             1.5174E-01
 GRADIENT:  -6.5361E-01  5.8569E-01  5.7755E-01  9.1530E-01 -1.1513E+00  5.9250E-02  7.4594E-02  1.4960E-02 -1.1854E-01  2.4340E-01
             6.4678E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1709.50677652488        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      758
 NPARAMETR:  9.9083E-01  5.5153E-01  1.0939E+00  1.3057E+00  8.5159E-01  9.9300E-01  7.9878E-01  2.8632E-02  8.7688E-01  7.9411E-01
             1.0523E+00
 PARAMETER:  9.0789E-02 -4.9506E-01  1.8978E-01  3.6678E-01 -6.0652E-02  9.2972E-02 -1.2467E-01 -3.4532E+00 -3.1385E-02 -1.3054E-01
             1.5098E-01
 GRADIENT:  -5.9608E-02 -7.4625E-02 -2.8755E-02 -1.0303E-01  2.7687E-02 -4.0705E-03 -4.9781E-04  4.7116E-03  3.7151E-03  6.7561E-03
             1.5068E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1709.50688038874        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      939
 NPARAMETR:  9.9084E-01  5.4657E-01  1.0940E+00  1.3086E+00  8.5012E-01  9.9289E-01  8.0229E-01  2.4486E-02  8.7517E-01  7.9403E-01
             1.0522E+00
 PARAMETER:  9.0793E-02 -5.0409E-01  1.8980E-01  3.6895E-01 -6.2374E-02  9.2865E-02 -1.2028E-01 -3.6096E+00 -3.3339E-02 -1.3063E-01
             1.5089E-01
 GRADIENT:   9.8347E-02 -2.0567E-01 -1.6268E-01 -3.3606E-01  3.0242E-01 -1.9692E-02  3.6020E-03  3.4484E-03  5.0557E-02 -8.8407E-04
             2.5790E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1709.51028925569        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1117
 NPARAMETR:  9.9114E-01  5.7922E-01  1.0946E+00  1.2893E+00  8.6081E-01  9.9340E-01  7.7956E-01  1.0000E-02  8.8730E-01  7.9545E-01
             1.0530E+00
 PARAMETER:  9.1105E-02 -4.4607E-01  1.9036E-01  3.5408E-01 -4.9884E-02  9.3381E-02 -1.4903E-01 -5.2303E+00 -1.9570E-02 -1.2885E-01
             1.5165E-01
 GRADIENT:  -1.3015E-01  3.2686E-02  2.5804E-02  3.0282E-02 -8.6795E-02  1.8628E-02 -8.1717E-04  0.0000E+00 -2.1652E-02  1.7876E-02
             7.1266E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1709.51030410109        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1294
 NPARAMETR:  9.9121E-01  5.7985E-01  1.0948E+00  1.2889E+00  8.6112E-01  9.9336E-01  7.7881E-01  1.0000E-02  8.8758E-01  7.9552E-01
             1.0530E+00
 PARAMETER:  9.1175E-02 -4.4515E-01  1.9057E-01  3.5380E-01 -4.9485E-02  9.3342E-02 -1.4984E-01 -5.2287E+00 -1.9232E-02 -1.2876E-01
             1.5168E-01
 GRADIENT:   3.0353E-03 -3.4662E-02 -3.1970E-02 -4.7089E-02  6.2694E-02  2.7285E-04  4.4305E-04  0.0000E+00  6.3285E-03  5.9751E-04
             7.6874E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1294
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.8409E-04 -1.1615E-02 -1.9645E-04 -1.9627E-03 -2.0688E-02
 SE:             2.9813E-02  8.7902E-03  2.0216E-04  2.8224E-02  2.4481E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9507E-01  1.8639E-01  3.3119E-01  9.4456E-01  3.9807E-01

 ETASHRINKSD(%)  1.2416E-01  7.0552E+01  9.9323E+01  5.4462E+00  1.7985E+01
 ETASHRINKVR(%)  2.4817E-01  9.1328E+01  9.9995E+01  1.0596E+01  3.2735E+01
 EBVSHRINKSD(%)  4.6689E-01  7.1132E+01  9.9305E+01  5.3261E+00  1.6865E+01
 EBVSHRINKVR(%)  9.3159E-01  9.1666E+01  9.9995E+01  1.0368E+01  3.0885E+01
 RELATIVEINF(%)  9.7840E+01  2.7217E-01  6.6474E-04  4.6977E+00  3.6294E+00
 EPSSHRINKSD(%)  4.1139E+01
 EPSSHRINKVR(%)  6.5354E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1709.5103041010884     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -974.35947753735024     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.93
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1709.510       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  5.80E-01  1.09E+00  1.29E+00  8.61E-01  9.93E-01  7.79E-01  1.00E-02  8.88E-01  7.96E-01  1.05E+00
 


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
+        1.14E+03
 
 TH 2
+       -1.97E+01  4.68E+02
 
 TH 3
+        1.08E+01  2.30E+02  4.55E+02
 
 TH 4
+       -1.08E+01  4.51E+02 -4.73E+01  7.61E+02
 
 TH 5
+       -7.79E-01 -5.34E+02 -7.88E+02 -4.88E+00  1.64E+03
 
 TH 6
+        6.05E-01 -4.32E+00  2.58E+00 -2.31E+00 -3.86E+00  1.98E+02
 
 TH 7
+       -6.79E-01 -3.50E+00  5.04E+00 -5.12E+00 -2.06E+00 -2.37E+00  5.48E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.45E+00 -6.10E+01  7.06E+00  1.13E+01 -6.73E+00 -1.94E+00  1.60E+01  0.00E+00  2.05E+02
 
 TH10
+        6.99E-01  2.10E+01 -1.82E+01 -5.32E+00 -8.20E+01  5.90E-02  9.21E+00  0.00E+00  5.09E-01  1.51E+02
 
 TH11
+       -5.59E+00 -1.86E+01 -3.60E+01 -8.60E+00  1.20E+01  4.05E-01  1.45E+00  0.00E+00  8.65E+00  4.09E+01  2.09E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.153
Stop Time:
Sat Sep 18 11:44:10 CDT 2021

Thu Sep 30 07:43:49 CDT 2021
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
$DATA ../../../../data/spa2/TD1/dat87.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2167.64961662585        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1890E+02 -3.7936E+01  2.7028E+01 -2.6625E+01  1.2725E+02  4.2865E+01  4.2244E+00 -4.6226E+02 -1.1701E+02  7.6200E+00
            -8.4287E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2354.10354270820        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0113E+00  9.8975E-01  1.1654E+00  1.0166E+00  1.0160E+00  1.0350E+00  9.4867E-01  1.9181E+00  8.2356E-01  9.2119E-01
             9.5675E-01
 PARAMETER:  1.1125E-01  8.9700E-02  2.5304E-01  1.1651E-01  1.1584E-01  1.3437E-01  4.7305E-02  7.5133E-01 -9.4122E-02  1.7908E-02
             5.5789E-02
 GRADIENT:   5.5203E+02 -1.4610E+01  5.9824E-01 -2.5642E+01  5.3171E+01  8.9187E+01 -1.1928E+01 -1.3067E+02 -5.5609E+01 -8.3230E+00
            -1.0154E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2373.32827558514        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      243
 NPARAMETR:  1.0112E+00  9.6256E-01  1.9249E+00  1.1276E+00  1.2175E+00  1.0837E+00  1.3356E+00  2.1613E+00  8.3588E-01  1.0778E+00
             1.0530E+00
 PARAMETER:  1.1113E-01  6.1839E-02  7.5485E-01  2.2007E-01  2.9676E-01  1.8037E-01  3.8938E-01  8.7070E-01 -7.9274E-02  1.7491E-01
             1.5166E-01
 GRADIENT:  -3.0811E+00  1.2410E+00  2.1342E+01  6.0570E+00  6.4291E+01  1.7082E+01  8.5963E+00 -1.6156E+02 -1.1469E+01 -4.7110E+00
             2.7623E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2409.16278555704        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      420
 NPARAMETR:  1.0135E+00  9.1631E-01  1.7655E+00  1.1419E+00  1.1233E+00  1.0445E+00  1.2358E+00  3.2057E+00  8.3368E-01  9.6882E-01
             1.0276E+00
 PARAMETER:  1.1345E-01  1.2599E-02  6.6845E-01  2.3272E-01  2.1626E-01  1.4355E-01  3.1172E-01  1.2649E+00 -8.1907E-02  6.8321E-02
             1.2720E-01
 GRADIENT:   1.6388E+00 -1.4806E+01 -2.6003E+00 -4.5067E+00  1.6257E+01  2.9693E+00 -3.6088E+00 -1.5081E+01 -9.4798E+00 -7.2379E+00
            -8.8928E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2410.40385149184        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      608             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0136E+00  9.8284E-01  1.7853E+00  1.1072E+00  1.1570E+00  1.0386E+00  1.1928E+00  3.3161E+00  8.7507E-01  1.0159E+00
             1.0407E+00
 PARAMETER:  1.1356E-01  8.2695E-02  6.7961E-01  2.0185E-01  2.4585E-01  1.3791E-01  2.7630E-01  1.2988E+00 -3.3456E-02  1.1576E-01
             1.3988E-01
 GRADIENT:   4.8744E+02  4.8712E+01  1.9746E+01  1.8726E+02  5.8384E+01  8.2295E+01  3.9994E+01  4.0124E+01  3.7241E+00 -1.1087E+00
             1.6431E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2410.49234705137        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      796
 NPARAMETR:  1.0134E+00  9.8522E-01  1.7853E+00  1.1073E+00  1.1536E+00  1.0384E+00  1.1886E+00  3.3243E+00  8.7976E-01  1.0217E+00
             1.0411E+00
 PARAMETER:  1.1333E-01  8.5114E-02  6.7960E-01  2.0189E-01  2.4292E-01  1.3770E-01  2.7276E-01  1.3013E+00 -2.8105E-02  1.2147E-01
             1.4024E-01
 GRADIENT:   8.2313E-01 -4.2233E+00 -1.4422E+00  5.4758E+00  6.9753E+00  5.3392E-01 -5.9589E-02 -6.8615E+00 -2.9352E+00 -2.2587E+00
             1.9182E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2410.52354851431        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      984
 NPARAMETR:  1.0132E+00  9.8663E-01  1.7853E+00  1.1071E+00  1.1516E+00  1.0383E+00  1.1848E+00  3.3247E+00  8.8281E-01  1.0244E+00
             1.0414E+00
 PARAMETER:  1.1309E-01  8.6541E-02  6.7960E-01  2.0173E-01  2.4114E-01  1.3754E-01  2.6954E-01  1.3014E+00 -2.4645E-02  1.2406E-01
             1.4056E-01
 GRADIENT:   3.2515E-01 -3.1787E+00 -1.0402E+00  6.3608E+00  4.3041E+00  4.7384E-01 -5.6407E-02 -6.5286E+00 -2.5695E+00 -1.4894E+00
             4.3741E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2410.54849417413        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1141
 NPARAMETR:  1.0138E+00  9.8613E-01  1.7853E+00  1.1060E+00  1.1530E+00  1.0375E+00  1.1726E+00  3.3255E+00  8.9910E-01  1.0250E+00
             1.0411E+00
 PARAMETER:  1.1366E-01  8.6037E-02  6.7959E-01  2.0071E-01  2.4235E-01  1.3679E-01  2.5921E-01  1.3016E+00 -6.3566E-03  1.2469E-01
             1.4030E-01
 GRADIENT:   1.5570E+00 -5.9180E+00 -1.1629E+00  4.7924E+00  6.0914E+00  1.7501E-01  6.4026E-02 -6.3201E+00  6.6973E-02 -1.2582E+00
             1.5894E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2410.54849417413        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  1.0138E+00  9.8613E-01  1.7853E+00  1.1060E+00  1.1530E+00  1.0375E+00  1.1726E+00  3.3255E+00  8.9910E-01  1.0250E+00
             1.0411E+00
 PARAMETER:  1.1366E-01  8.6037E-02  6.7959E-01  2.0071E-01  2.4235E-01  1.3679E-01  2.5921E-01  1.3016E+00 -6.3566E-03  1.2469E-01
             1.4030E-01
 GRADIENT:   1.5570E+00 -5.9180E+00 -1.1629E+00  4.7924E+00  6.0914E+00  1.7501E-01  6.4026E-02 -6.3201E+00  6.6973E-02 -1.2582E+00
             1.5894E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1163
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1112E-03 -1.4046E-02 -4.4665E-02  3.5245E-03 -6.7238E-02
 SE:             2.9914E-02  1.9619E-02  2.4288E-02  2.4676E-02  2.0504E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7037E-01  4.7403E-01  6.5923E-02  8.8642E-01  1.0408E-03

 ETASHRINKSD(%)  1.0000E-10  3.4275E+01  1.8631E+01  1.7334E+01  3.1309E+01
 ETASHRINKVR(%)  1.0000E-10  5.6803E+01  3.3790E+01  3.1663E+01  5.2816E+01
 EBVSHRINKSD(%)  3.3404E-01  3.4412E+01  1.9486E+01  1.9329E+01  3.0583E+01
 EBVSHRINKVR(%)  6.6696E-01  5.6982E+01  3.5175E+01  3.4922E+01  5.1813E+01
 RELATIVEINF(%)  9.9310E+01  1.0847E+01  4.0155E+01  1.7715E+01  2.8393E+01
 EPSSHRINKSD(%)  3.2307E+01
 EPSSHRINKVR(%)  5.4176E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2410.5484941741333     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1307.8222543285262     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.85
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2410.548       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  9.86E-01  1.79E+00  1.11E+00  1.15E+00  1.04E+00  1.17E+00  3.33E+00  8.99E-01  1.02E+00  1.04E+00
 


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
+        3.04E+06
 
 TH 2
+        3.14E+03  4.15E+06
 
 TH 3
+        5.71E-01  1.85E+01  5.49E+04
 
 TH 4
+        3.79E+03 -1.23E+03 -2.36E+01  1.64E+06
 
 TH 5
+        1.25E+06 -1.47E+06 -4.75E+01 -1.30E+06  1.04E+06
 
 TH 6
+       -9.33E+01  1.08E+02  1.28E-01  4.79E+01 -3.96E+01  1.83E+02
 
 TH 7
+       -4.50E+02  5.55E+02 -1.28E+00  2.23E+02 -1.99E+02 -4.95E-02  3.80E+01
 
 TH 8
+       -8.32E+04  9.72E+04 -5.94E+00  4.31E+04 -3.43E+04 -3.02E-02  1.20E-01  4.33E+03
 
 TH 9
+       -1.13E+03  1.31E+03  3.00E-01  6.00E+02 -4.58E+02 -1.78E-01 -1.48E+06  2.75E+00  5.00E+06
 
 TH10
+        2.74E+06 -6.41E+06 -3.22E+00 -1.43E+06  1.13E+06 -8.46E+01 -4.11E+02 -7.50E+04 -1.01E+03  4.95E+06
 
 TH11
+       -1.01E+02  7.46E+01  5.22E-01  1.25E+06 -9.90E+05  1.59E+00  4.11E+00  6.56E+04 -1.37E+00 -2.16E+06  1.89E+06
 
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
 #CPUT: Total CPU Time in Seconds,       35.098
Stop Time:
Thu Sep 30 07:44:26 CDT 2021

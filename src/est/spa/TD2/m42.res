Sat Sep 18 14:40:20 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat42.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1705.08104283113        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -6.3577E+01  1.5142E+01 -2.6112E+01  5.4874E+01  6.2407E+01 -2.6295E-01  4.0905E+00  2.1039E+00 -3.6765E+00  4.5992E+00
            -3.8593E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1710.92204672629        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0481E+00  8.9412E-01  9.3960E-01  1.0375E+00  8.8010E-01  9.9992E-01  9.0675E-01  9.8107E-01  1.0721E+00  8.2058E-01
             1.0440E+00
 PARAMETER:  1.4694E-01 -1.1910E-02  3.7701E-02  1.3678E-01 -2.7719E-02  9.9922E-02  2.1151E-03  8.0884E-02  1.6965E-01 -9.7742E-02
             1.4305E-01
 GRADIENT:  -1.0738E+01 -1.6372E+00 -3.6164E+00  9.1679E+00  7.5938E+00 -8.8397E-02  2.6873E+00  1.9935E+00  1.6383E+01 -4.1232E+00
             1.3865E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1711.55481698232        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  1.0520E+00  7.9893E-01  9.2501E-01  1.1059E+00  8.3682E-01  1.0005E+00  8.5726E-01  8.6462E-01  1.0109E+00  8.3364E-01
             1.0255E+00
 PARAMETER:  1.5074E-01 -1.2448E-01  2.2047E-02  2.0066E-01 -7.8140E-02  1.0050E-01 -5.4011E-02 -4.5463E-02  1.1080E-01 -8.1950E-02
             1.2519E-01
 GRADIENT:  -1.5561E+00  8.4067E+00 -6.5701E+00  2.9883E+01  4.8773E+00  2.9656E-01 -5.8296E-01  4.1022E-01  9.1819E+00 -3.3125E+00
             6.4189E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1712.42710275170        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  1.0536E+00  8.3764E-01  8.0265E-01  1.0535E+00  7.9366E-01  1.0017E+00  1.1253E+00  6.2751E-01  9.3257E-01  8.1530E-01
             1.0052E+00
 PARAMETER:  1.5220E-01 -7.7165E-02 -1.1983E-01  1.5207E-01 -1.3110E-01  1.0170E-01  2.1806E-01 -3.6599E-01  3.0186E-02 -1.0420E-01
             1.0521E-01
 GRADIENT:  -1.0229E+00 -2.3562E+00 -1.1313E+00 -4.9284E+00  7.7330E-01  2.4048E-01  6.1754E-01  9.0828E-01 -7.1948E-01  1.3176E+00
            -6.2068E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1712.66954648854        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  1.0569E+00  9.5646E-01  6.7700E-01  9.7489E-01  7.7568E-01  1.0025E+00  1.0488E+00  3.2073E-01  9.7322E-01  7.9254E-01
             1.0071E+00
 PARAMETER:  1.5539E-01  5.5484E-02 -2.9008E-01  7.4569E-02 -1.5401E-01  1.0245E-01  1.4762E-01 -1.0372E+00  7.2859E-02 -1.3251E-01
             1.0711E-01
 GRADIENT:   1.8814E+00  1.1675E+00 -1.6709E+00  5.0294E+00  1.0330E+00 -2.3917E-01 -1.6930E-01  2.1629E-01  7.5906E-01  3.8940E-01
             2.4765E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1712.75968962987        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  1.0565E+00  1.0434E+00  6.2874E-01  9.1460E-01  7.8805E-01  1.0034E+00  9.9437E-01  1.5939E-01  1.0076E+00  7.8429E-01
             1.0063E+00
 PARAMETER:  1.5494E-01  1.4251E-01 -3.6404E-01  1.0735E-02 -1.3820E-01  1.0341E-01  9.4359E-02 -1.7364E+00  1.0752E-01 -1.4297E-01
             1.0630E-01
 GRADIENT:  -1.1955E-01 -1.1137E+00 -9.4399E-01 -9.7307E-01  1.3532E+00 -8.6768E-02  6.8659E-02  7.7245E-02 -5.4107E-03 -6.4584E-02
             9.8069E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1712.79613266376        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1055
 NPARAMETR:  1.0564E+00  1.0182E+00  6.2954E-01  9.3005E-01  7.7657E-01  1.0035E+00  1.0145E+00  3.8563E-02  9.9309E-01  7.8372E-01
             1.0063E+00
 PARAMETER:  1.5489E-01  1.1801E-01 -3.6277E-01  2.7481E-02 -1.5286E-01  1.0354E-01  1.1443E-01 -3.1555E+00  9.3065E-02 -1.4370E-01
             1.0631E-01
 GRADIENT:  -2.4331E-01 -2.1315E-01 -4.8115E-01  8.0990E-02  1.1395E-01 -2.8661E-02  3.8344E-02  4.0935E-03 -5.2154E-02  1.4453E-01
             1.6663E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1712.79886011735        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1230
 NPARAMETR:  1.0565E+00  1.0234E+00  6.3082E-01  9.2719E-01  7.7983E-01  1.0036E+00  1.0101E+00  1.0000E-02  9.9665E-01  7.8604E-01
             1.0061E+00
 PARAMETER:  1.5498E-01  1.2316E-01 -3.6074E-01  2.4405E-02 -1.4868E-01  1.0360E-01  1.1010E-01 -4.6301E+00  9.6640E-02 -1.4075E-01
             1.0611E-01
 GRADIENT:   1.0104E-02  6.9406E-03 -1.1025E-03  1.2279E-02  2.5344E-03 -1.2801E-03 -4.7860E-03  0.0000E+00 -5.4096E-03 -3.2827E-03
            -4.7873E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1712.79886011735        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  1.0565E+00  1.0234E+00  6.3082E-01  9.2719E-01  7.7983E-01  1.0036E+00  1.0101E+00  1.0000E-02  9.9665E-01  7.8604E-01
             1.0061E+00
 PARAMETER:  1.5498E-01  1.2316E-01 -3.6074E-01  2.4405E-02 -1.4868E-01  1.0360E-01  1.1010E-01 -4.6301E+00  9.6640E-02 -1.4075E-01
             1.0611E-01
 GRADIENT:   1.0104E-02  6.9406E-03 -1.1025E-03  1.2279E-02  2.5344E-03 -1.2801E-03 -4.7860E-03  0.0000E+00 -5.4096E-03 -3.2827E-03
            -4.7873E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1252
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.9109E-06 -9.8622E-03 -4.0866E-04  3.8775E-03 -1.5503E-02
 SE:             2.9851E-02  2.0758E-02  1.8471E-04  2.5630E-02  2.2855E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9984E-01  6.3472E-01  2.6936E-02  8.7975E-01  4.9756E-01

 ETASHRINKSD(%)  1.0000E-10  3.0457E+01  9.9381E+01  1.4135E+01  2.3433E+01
 ETASHRINKVR(%)  1.0000E-10  5.1638E+01  9.9996E+01  2.6272E+01  4.1376E+01
 EBVSHRINKSD(%)  4.2870E-01  3.0178E+01  9.9439E+01  1.4277E+01  2.3158E+01
 EBVSHRINKVR(%)  8.5557E-01  5.1249E+01  9.9997E+01  2.6516E+01  4.0953E+01
 RELATIVEINF(%)  9.8973E+01  2.4940E+00  3.1130E-04  5.2449E+00  4.5939E+00
 EPSSHRINKSD(%)  4.4467E+01
 EPSSHRINKVR(%)  6.9161E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1712.7988601173502     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -977.64803355361198     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.38
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1712.799       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.02E+00  6.31E-01  9.27E-01  7.80E-01  1.00E+00  1.01E+00  1.00E-02  9.97E-01  7.86E-01  1.01E+00
 


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
+        9.80E+02
 
 TH 2
+       -6.60E+00  5.29E+02
 
 TH 3
+        1.47E+01  3.63E+02  9.34E+02
 
 TH 4
+       -1.34E+01  3.55E+02 -3.71E+02  9.63E+02
 
 TH 5
+       -2.04E+00 -5.90E+02 -1.11E+03  3.76E+02  1.73E+03
 
 TH 6
+       -2.20E-01 -5.88E-01  2.79E+00 -2.63E+00  1.42E-01  1.96E+02
 
 TH 7
+        9.07E-01  2.40E+01 -1.02E+01 -6.52E+00 -1.66E+01  9.62E-01  5.01E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.27E+00 -2.75E+01 -3.23E+01  3.59E+01 -2.79E+00  2.06E-01  2.19E+01  0.00E+00  1.12E+02
 
 TH10
+       -4.90E-01 -9.29E+00 -8.59E+01 -1.95E+01 -5.09E+01 -6.61E-01  2.08E+01  0.00E+00  6.01E+00  1.24E+02
 
 TH11
+       -5.97E+00 -1.44E+01 -3.84E+01 -6.05E+00  9.03E+00  1.22E+00  7.70E+00  0.00E+00  7.99E+00  2.40E+01  2.09E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.196
Stop Time:
Sat Sep 18 14:40:43 CDT 2021

Sat Sep 25 11:50:44 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat67.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1619.36099177157        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.4878E+01 -1.0406E+02 -6.0664E+00 -1.5397E+02  4.6915E+00 -4.6126E+00 -1.5394E+01  6.6883E+00 -2.2679E+01 -5.0348E+00
            -4.1808E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1631.97278914936        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.8043E-01  1.1140E+00  1.0873E+00  1.0391E+00  1.0772E+00  1.0145E+00  1.0996E+00  9.3326E-01  1.0482E+00  1.0600E+00
             1.1086E+00
 PARAMETER:  8.0237E-02  2.0796E-01  1.8372E-01  1.3836E-01  1.7440E-01  1.1440E-01  1.9490E-01  3.0928E-02  1.4710E-01  1.5823E-01
             2.0310E-01
 GRADIENT:   1.2099E+01  6.4120E+00  1.0073E+01 -2.3409E+00 -8.0498E+00  1.2745E+00  6.2004E-01 -6.6140E-01  1.8033E+00 -5.2698E+00
             3.3331E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1632.30292034441        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      208
 NPARAMETR:  9.8621E-01  1.1133E+00  1.0734E+00  1.0495E+00  1.0832E+00  1.0322E+00  1.1123E+00  8.9554E-01  1.0334E+00  1.1139E+00
             1.1026E+00
 PARAMETER:  8.6114E-02  2.0729E-01  1.7079E-01  1.4835E-01  1.7991E-01  1.3174E-01  2.0645E-01 -1.0327E-02  1.3287E-01  2.0788E-01
             1.9767E-01
 GRADIENT:  -7.8010E+00  3.8093E+00  1.7637E+00  4.7341E+00 -1.0314E+00  3.3159E+00  1.9456E-01 -8.2882E-03 -3.0175E-01  7.9640E-01
             2.0024E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1632.63852229287        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      385
 NPARAMETR:  9.9362E-01  1.2710E+00  8.6192E-01  9.3804E-01  1.0495E+00  1.0199E+00  1.0223E+00  6.0917E-01  1.1060E+00  1.0535E+00
             1.0913E+00
 PARAMETER:  9.3602E-02  3.3979E-01 -4.8589E-02  3.6034E-02  1.4830E-01  1.1974E-01  1.2208E-01 -3.9566E-01  2.0076E-01  1.5215E-01
             1.8736E-01
 GRADIENT:   4.1051E+00  3.4200E+00  1.4165E+00  3.6039E+00 -4.0440E+00 -2.1283E+00 -1.4490E+00  2.8489E-01 -7.3047E-01 -7.7355E-01
            -2.0279E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1633.20153687237        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      563
 NPARAMETR:  9.9264E-01  1.6261E+00  6.4490E-01  7.0176E-01  1.1433E+00  1.0275E+00  8.7964E-01  2.2263E-01  1.3551E+00  1.0844E+00
             1.0965E+00
 PARAMETER:  9.2611E-02  5.8619E-01 -3.3866E-01 -2.5416E-01  2.3395E-01  1.2711E-01 -2.8239E-02 -1.4022E+00  4.0387E-01  1.8106E-01
             1.9216E-01
 GRADIENT:  -9.7929E-01  4.4893E+00 -1.4151E+00  2.8203E+00 -2.5800E+00  6.5596E-02  8.8424E-01  1.6616E-01  2.8527E-01  2.7163E-01
             2.4858E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1633.24904778879        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      738
 NPARAMETR:  9.9319E-01  1.6677E+00  6.3634E-01  6.7038E-01  1.1741E+00  1.0269E+00  8.5930E-01  1.8506E-01  1.4041E+00  1.1065E+00
             1.0973E+00
 PARAMETER:  9.3169E-02  6.1144E-01 -3.5203E-01 -2.9992E-01  2.6050E-01  1.2650E-01 -5.1635E-02 -1.5871E+00  4.3941E-01  2.0122E-01
             1.9283E-01
 GRADIENT:   2.9245E-01 -1.5516E+00 -4.7383E-01 -7.6455E-01  3.4002E-01 -1.7282E-01  9.0327E-02  9.3374E-02  1.3495E-01  1.1013E-01
             3.8272E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1633.29575072964        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      916
 NPARAMETR:  9.9307E-01  1.6335E+00  6.4828E-01  6.9399E-01  1.1566E+00  1.0271E+00  8.7325E-01  4.4398E-02  1.3689E+00  1.0968E+00
             1.0967E+00
 PARAMETER:  9.3042E-02  5.9074E-01 -3.3344E-01 -2.6530E-01  2.4551E-01  1.2671E-01 -3.5535E-02 -3.0146E+00  4.1397E-01  1.9244E-01
             1.9232E-01
 GRADIENT:   7.4515E-02  5.7532E-01  1.6340E-01  2.9453E-01  9.9660E-02 -7.4196E-02  3.0924E-02  4.7863E-03 -1.3733E-01 -7.5865E-02
            -9.4405E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1633.29879136271        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1091
 NPARAMETR:  9.9304E-01  1.6435E+00  6.4309E-01  6.8699E-01  1.1604E+00  1.0273E+00  8.6897E-01  1.0000E-02  1.3794E+00  1.0987E+00
             1.0969E+00
 PARAMETER:  9.3020E-02  5.9680E-01 -3.4148E-01 -2.7544E-01  2.4875E-01  1.2692E-01 -4.0451E-02 -4.6130E+00  4.2166E-01  1.9413E-01
             1.9251E-01
 GRADIENT:  -1.1106E-03 -2.4412E-02 -5.7265E-03 -9.0444E-03 -4.8429E-05  1.3824E-03  1.1867E-03  0.0000E+00  3.6653E-03  1.9058E-03
             2.6995E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1633.29879136271        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1113
 NPARAMETR:  9.9304E-01  1.6435E+00  6.4309E-01  6.8699E-01  1.1604E+00  1.0273E+00  8.6897E-01  1.0000E-02  1.3794E+00  1.0987E+00
             1.0969E+00
 PARAMETER:  9.3020E-02  5.9680E-01 -3.4148E-01 -2.7544E-01  2.4875E-01  1.2692E-01 -4.0451E-02 -4.6130E+00  4.2166E-01  1.9413E-01
             1.9251E-01
 GRADIENT:  -1.1106E-03 -2.4412E-02 -5.7265E-03 -9.0444E-03 -4.8429E-05  1.3824E-03  1.1867E-03  0.0000E+00  3.6653E-03  1.9058E-03
             2.6995E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1113
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.6344E-05 -2.9701E-02 -3.3142E-04  2.1330E-02 -3.8490E-02
 SE:             2.9812E-02  2.2541E-02  1.2228E-04  2.3477E-02  2.2277E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9796E-01  1.8761E-01  6.7227E-03  3.6358E-01  8.4026E-02

 ETASHRINKSD(%)  1.2487E-01  2.4486E+01  9.9590E+01  2.1350E+01  2.5369E+01
 ETASHRINKVR(%)  2.4959E-01  4.2976E+01  9.9998E+01  3.8142E+01  4.4302E+01
 EBVSHRINKSD(%)  4.9181E-01  2.3363E+01  9.9650E+01  2.2802E+01  2.3845E+01
 EBVSHRINKVR(%)  9.8120E-01  4.1268E+01  9.9999E+01  4.0405E+01  4.2005E+01
 RELATIVEINF(%)  9.8937E+01  3.9428E+00  1.8077E-04  4.4740E+00  1.1829E+01
 EPSSHRINKSD(%)  4.3595E+01
 EPSSHRINKVR(%)  6.8185E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1633.2987913627103     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -898.14796479897211     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1633.299       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.93E-01  1.64E+00  6.43E-01  6.87E-01  1.16E+00  1.03E+00  8.69E-01  1.00E-02  1.38E+00  1.10E+00  1.10E+00
 


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
+        1.06E+03
 
 TH 2
+       -5.13E+00  3.40E+02
 
 TH 3
+        7.78E+00  1.57E+02  3.84E+02
 
 TH 4
+       -1.48E+01  2.70E+02 -2.22E+02  7.47E+02
 
 TH 5
+       -3.61E+00 -1.61E+02 -2.63E+02  1.81E+02  3.97E+02
 
 TH 6
+       -1.09E+00 -8.00E-01  1.00E+00 -5.23E+00 -2.34E-01  1.86E+02
 
 TH 7
+        4.89E+00  8.93E+00 -1.27E+01 -8.62E+00 -9.56E+00  3.47E+00  1.04E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.33E-01 -1.88E+01 -2.97E+01  4.57E+01  3.39E-01 -2.86E-01  1.82E+01  0.00E+00  4.48E+01
 
 TH10
+       -2.34E+00 -1.09E+01 -3.07E+01 -5.90E+00 -4.91E+01  1.01E+00  8.89E+00  0.00E+00  5.17E+00  6.58E+01
 
 TH11
+       -8.08E+00 -1.80E+01 -3.03E+01  7.00E-01 -1.83E+00  2.31E+00  8.87E+00  0.00E+00  5.96E+00  1.22E+01  1.82E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.820
Stop Time:
Sat Sep 25 11:51:05 CDT 2021

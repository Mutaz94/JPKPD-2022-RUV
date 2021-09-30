Thu Sep 30 00:06:42 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat26.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -209.670999314240        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3924E+02  7.6273E+00  1.8392E+02 -4.7172E+01  1.1089E+02  6.4485E+01 -4.3975E+01 -2.4601E+02 -7.8678E+01 -1.2964E+02
            -3.2410E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1454.49164392958        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0572E+00  1.0895E+00  8.9047E-01  1.1263E+00  9.5139E-01  7.6754E-01  9.6686E-01  1.1661E+00  8.9749E-01  1.0485E+00
             5.2882E+00
 PARAMETER:  1.5562E-01  1.8569E-01 -1.6007E-02  2.1891E-01  5.0171E-02 -1.6457E-01  6.6299E-02  2.5366E-01 -8.1506E-03  1.4731E-01
             1.7655E+00
 GRADIENT:   9.1887E+01  1.2038E+01 -1.3348E+01  3.6204E+01 -1.7267E+01 -3.7140E+00  1.2355E+01  9.6023E+00  2.3952E+01  2.2176E+01
             3.4160E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1529.79968709959        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.7855E-01  1.3316E+00  9.1751E-01  8.9504E-01  1.0776E+00  7.5117E-01  2.7910E-01  1.7185E+00  1.2600E+00  5.3033E-01
             4.0256E+00
 PARAMETER:  7.8319E-02  3.8636E-01  1.3905E-02 -1.0891E-02  1.7474E-01 -1.8612E-01 -1.1762E+00  6.4147E-01  3.3110E-01 -5.3425E-01
             1.4927E+00
 GRADIENT:  -7.6779E+01  2.4378E+01 -3.7592E+00  2.3087E+01 -3.4229E+01 -2.4389E+01 -1.0874E+00  1.2141E+01  1.8181E+01  3.3219E+00
             1.7312E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1557.99958267547        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  9.8052E-01  1.3864E+00  1.6802E+00  8.7241E-01  1.5944E+00  8.4594E-01  8.7981E-01  2.4992E+00  1.0160E+00  4.2481E-01
             3.2974E+00
 PARAMETER:  8.0333E-02  4.2669E-01  6.1894E-01 -3.6495E-02  5.6649E-01 -6.7312E-02 -2.8049E-02  1.0160E+00  1.1592E-01 -7.5612E-01
             1.2931E+00
 GRADIENT:   2.5457E-01  1.6526E+00 -7.8999E+00 -9.1703E-01  2.5990E+01  5.9541E+00  1.1373E+01 -1.2457E+01  5.9487E+00  9.4869E-01
             3.3224E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1566.80955145995        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.8157E-01  1.1574E+00  4.1514E+00  1.0486E+00  1.6882E+00  8.1839E-01  5.4969E-01  4.2165E+00  9.8631E-01  2.5418E-01
             3.3203E+00
 PARAMETER:  8.1400E-02  2.4617E-01  1.5234E+00  1.4746E-01  6.2369E-01 -1.0042E-01 -4.9840E-01  1.5390E+00  8.6211E-02 -1.2697E+00
             1.3001E+00
 GRADIENT:  -1.5269E+00  9.9614E+00 -2.1868E+00  1.7883E+01  8.4128E+00 -8.6280E-01 -6.8442E-01  1.0898E+00  3.4129E+00  2.6016E-01
             3.1900E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1567.84676290990        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  9.8177E-01  1.1413E+00  4.7842E+00  1.0394E+00  1.6624E+00  8.2298E-01  7.1780E-01  4.5657E+00  8.9879E-01  1.4166E-01
             3.3127E+00
 PARAMETER:  8.1602E-02  2.3220E-01  1.6653E+00  1.3869E-01  6.0823E-01 -9.4822E-02 -2.3156E-01  1.6186E+00 -6.7025E-03 -1.8544E+00
             1.2978E+00
 GRADIENT:   3.0803E+00 -1.3398E+01 -2.2229E+00 -1.6029E+01  3.2943E+00  1.5230E+00 -1.8174E-01  4.3625E+00  1.0578E-01  7.6544E-02
            -2.3654E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1569.08577063949        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      492
 NPARAMETR:  9.9348E-01  1.1245E+00  5.8022E+00  1.0772E+00  1.6788E+00  8.2618E-01  9.2807E-01  5.1363E+00  8.0226E-01  7.8360E-02
             3.3658E+00
 PARAMETER:  9.3455E-02  2.1731E-01  1.8582E+00  1.7440E-01  6.1806E-01 -9.0947E-02  2.5350E-02  1.7363E+00 -1.2033E-01 -2.4464E+00
             1.3137E+00
 GRADIENT:  -2.7398E+00  1.2361E+00 -1.4583E+00  1.0990E+00  4.7450E-01  1.0206E+00  3.0051E+00 -3.6672E+00  6.2814E-01  2.5464E-02
             3.1251E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1570.73520705576        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  9.9447E-01  1.0042E+00  1.2308E+01  1.1635E+00  1.7015E+00  8.2040E-01  1.0400E+00  7.3652E+00  7.0542E-01  3.0567E-02
             3.3599E+00
 PARAMETER:  9.4451E-02  1.0418E-01  2.6102E+00  2.5142E-01  6.3150E-01 -9.7962E-02  1.3924E-01  2.0968E+00 -2.4896E-01 -3.3878E+00
             1.3119E+00
 GRADIENT:  -2.2286E+00 -4.3142E+00 -4.9843E+00  6.4025E+00 -3.8545E+00 -2.7140E-01 -2.7892E+00  9.9487E+00  3.0208E+00  4.8833E-03
            -9.6277E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1570.87224721788        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      836
 NPARAMETR:  9.9535E-01  9.8825E-01  1.4232E+01  1.1767E+00  1.7002E+00  8.2126E-01  1.0117E+00  7.5430E+00  7.3237E-01  2.6232E-02
             3.3549E+00
 PARAMETER:  9.5340E-02  8.8161E-02  2.7550E+00  2.6276E-01  6.3065E-01 -9.6921E-02  1.1162E-01  2.1209E+00 -2.1144E-01 -3.5057E+00
             1.3102E+00
 GRADIENT:  -1.7319E-01 -2.7391E-01 -1.2968E+01  2.1465E+00 -6.3284E+01 -2.1990E-02  4.0291E-02  1.4157E+01  3.6516E-01  4.1893E-03
            -3.0543E+01
 NUMSIGDIG:         3.2         2.3         2.4         2.7         2.3         3.2         2.8         2.5         2.5         0.6
                    2.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      836
 NO. OF SIG. DIGITS IN FINAL EST.:  0.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6408E-03 -5.7498E-03 -1.7044E-02 -5.2896E-03 -1.6748E-04
 SE:             2.8590E-02  1.5918E-02  9.1499E-03  2.0470E-02  2.6939E-04
 N:                     100         100         100         100         100

 P VAL.:         9.5424E-01  7.1795E-01  6.2495E-02  7.9609E-01  5.3415E-01

 ETASHRINKSD(%)  4.2186E+00  4.6672E+01  6.9347E+01  3.1424E+01  9.9098E+01
 ETASHRINKVR(%)  8.2592E+00  7.1561E+01  9.0604E+01  5.2973E+01  9.9992E+01
 EBVSHRINKSD(%)  4.1394E+00  4.5646E+01  8.0552E+01  3.3000E+01  9.8977E+01
 EBVSHRINKVR(%)  8.1075E+00  7.0457E+01  9.6218E+01  5.5110E+01  9.9990E+01
 RELATIVEINF(%)  9.0587E+01  1.0855E+00  2.9789E+00  1.6593E+00  7.9614E-03
 EPSSHRINKSD(%)  1.7970E+01
 EPSSHRINKVR(%)  3.2712E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1570.8722472178831     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -651.93371401321042     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.21
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1570.872       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  9.88E-01  1.42E+01  1.18E+00  1.70E+00  8.21E-01  1.01E+00  7.54E+00  7.32E-01  2.72E-02  3.35E+00
 


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
+        1.53E+03
 
 TH 2
+       -7.01E+01  5.24E+02
 
 TH 3
+       -4.92E-01  1.01E-01  7.00E-01
 
 TH 4
+       -6.34E+01 -3.21E+04 -3.27E+00  1.13E+04
 
 TH 5
+       -7.95E+00 -4.44E+01  1.46E+00 -3.14E+03  9.62E+02
 
 TH 6
+        3.29E+00 -1.54E+01 -2.96E-02 -1.56E+01 -2.43E-01  2.44E+02
 
 TH 7
+       -3.20E+00  1.39E+02  1.24E-01 -2.85E+04 -8.83E+00 -2.54E+00  1.15E+02
 
 TH 8
+        8.04E-01 -4.98E-01  3.49E+00  6.32E+00 -3.25E+00  1.52E-01 -1.91E-01 -3.45E+00
 
 TH 9
+        8.57E+00 -3.77E+02 -3.38E+00  2.11E+04 -1.07E+02 -3.09E+00 -2.23E+02  7.06E+00  4.20E+04
 
 TH10
+       -5.53E-02  5.42E-01  7.33E-03 -3.62E+00  1.48E-01 -4.82E-02  4.11E-01 -1.64E-02 -1.09E+00  5.86E+00
 
 TH11
+       -2.42E+01 -1.04E+01  6.78E+00 -3.62E+01  1.02E+01  6.87E+00  4.16E+00  2.83E+01 -1.37E+01  1.69E-01  1.18E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.469
Stop Time:
Thu Sep 30 00:07:03 CDT 2021

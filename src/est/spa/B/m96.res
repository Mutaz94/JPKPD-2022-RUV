Sat Sep 18 08:50:19 CDT 2021
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
$DATA ../../../../data/spa/B/dat96.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1650.98660299086        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1424E+02 -1.3550E+01 -4.8929E+01  4.2385E+01  4.8808E+01 -1.1727E+01  1.1894E+01  1.3421E+01  2.3058E+01  2.0554E+01
            -3.0215E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1659.38709681913        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.7763E-01  1.0322E+00  1.1313E+00  9.8707E-01  1.0343E+00  1.0217E+00  9.0040E-01  8.9444E-01  8.9060E-01  8.6281E-01
             1.0913E+00
 PARAMETER:  7.7379E-02  1.3170E-01  2.2335E-01  8.6984E-02  1.3377E-01  1.2142E-01 -4.9175E-03 -1.1563E-02 -1.5854E-02 -4.7559E-02
             1.8734E-01
 GRADIENT:   5.6475E+01  2.0094E+01 -6.9369E+00  3.5903E+01  2.8943E+01  1.7725E-01  1.5207E-01  3.5253E+00 -6.9994E+00 -9.2089E+00
            -1.8083E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1661.15622939689        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.7322E-01  9.7823E-01  1.1622E+00  1.0138E+00  1.0404E+00  1.0577E+00  6.9808E-01  5.0406E-01  9.4966E-01  1.0656E+00
             1.0930E+00
 PARAMETER:  7.2854E-02  7.7989E-02  2.5027E-01  1.1366E-01  1.3960E-01  1.5608E-01 -2.5943E-01 -5.8505E-01  4.8352E-02  1.6355E-01
             1.8888E-01
 GRADIENT:   4.7018E+01  9.1426E+00 -1.1510E+01  2.6638E+01  1.7989E+01  1.5397E+01  7.6445E-02  1.5481E+00 -4.6032E+00  1.0406E+01
             5.7036E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1662.22500757785        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.5604E-01  1.0092E+00  9.8883E-01  9.7570E-01  9.6704E-01  1.0262E+00  8.8033E-01  3.4130E-01  9.3283E-01  9.2300E-01
             1.0763E+00
 PARAMETER:  5.5044E-02  1.0915E-01  8.8763E-02  7.5396E-02  6.6480E-02  1.2583E-01 -2.7454E-02 -9.7500E-01  3.0472E-02  1.9876E-02
             1.7357E-01
 GRADIENT:   9.7661E+00 -1.8605E+00 -5.7701E+00  2.6900E+00  7.7536E+00  2.2091E+00  9.3951E-01  8.8759E-01  7.1262E-01  2.7795E+00
            -3.4155E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1662.29287107159        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  9.5074E-01  9.6368E-01  9.7136E-01  1.0019E+00  9.3376E-01  1.0208E+00  9.2632E-01  2.0828E-01  9.0411E-01  8.9103E-01
             1.0784E+00
 PARAMETER:  4.9487E-02  6.3000E-02  7.0944E-02  1.0187E-01  3.1464E-02  1.2055E-01  2.3470E-02 -1.4689E+00 -8.0183E-04 -1.5378E-02
             1.7546E-01
 GRADIENT:  -2.3905E+00  3.6943E-01 -1.0364E-01  1.1427E-01 -9.2853E-01 -4.4015E-01  1.8118E-01  3.0290E-01  1.5690E-01  3.4223E-01
             7.2471E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1662.52302073882        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.5348E-01  9.4478E-01  9.9116E-01  1.0160E+00  9.3615E-01  1.0211E+00  9.2874E-01  5.4945E-02  8.9697E-01  8.9143E-01
             1.0773E+00
 PARAMETER:  5.2367E-02  4.3193E-02  9.1120E-02  1.1584E-01  3.4024E-02  1.2084E-01  2.6071E-02 -2.8014E+00 -8.7336E-03 -1.4926E-02
             1.7450E-01
 GRADIENT:   4.2493E+00  1.8218E+00  2.1086E+00  4.0331E+00  1.3761E+00 -6.4766E-02 -3.6645E-01  1.5417E-02 -4.4783E-01 -2.6019E+00
            -1.5397E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1662.84080907402        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      552
 NPARAMETR:  9.6830E-01  9.2675E-01  9.8110E-01  1.0287E+00  9.2386E-01  1.0350E+00  9.5309E-01  2.3215E-02  8.8738E-01  8.9409E-01
             1.0791E+00
 PARAMETER:  6.7789E-02  2.3928E-02  8.0914E-02  1.2825E-01  2.0807E-02  1.3443E-01  5.1950E-02 -3.6630E+00 -1.9482E-02 -1.1953E-02
             1.7614E-01
 GRADIENT:   2.3041E-01  4.8506E-02 -1.7593E-01  8.0269E-02  2.7859E-01  4.7685E-01 -1.9383E-01  3.3186E-03 -2.9721E-01 -5.1957E-02
             1.2613E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1662.84968930561        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      728
 NPARAMETR:  9.6776E-01  8.4981E-01  9.9183E-01  1.0764E+00  8.9698E-01  1.0336E+00  1.0376E+00  1.0000E-02  8.5123E-01  8.8311E-01
             1.0784E+00
 PARAMETER:  6.7230E-02 -6.2747E-02  9.1794E-02  1.7363E-01 -8.7234E-03  1.3301E-01  1.3687E-01 -4.5852E+00 -6.1073E-02 -2.4302E-02
             1.7544E-01
 GRADIENT:   3.9895E-02  1.4925E-01 -4.7611E-03  2.7828E-01 -5.9617E-03  5.1549E-02 -3.2472E-03  0.0000E+00 -1.5038E-02  3.6565E-03
            -1.3678E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1662.84973039763        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:      862
 NPARAMETR:  9.6778E-01  8.4531E-01  9.9218E-01  1.0789E+00  8.9544E-01  1.0335E+00  1.0425E+00  1.0000E-02  8.4937E-01  8.8243E-01
             1.0784E+00
 PARAMETER:  6.7177E-02 -6.7955E-02  9.2218E-02  1.7602E-01 -1.0470E-02  1.3285E-01  1.4180E-01 -4.6527E+00 -6.3348E-02 -2.5110E-02
             1.7543E-01
 GRADIENT:  -2.3054E-02  1.1105E-02  7.8797E-03  1.3491E-02 -1.0309E-02 -7.8432E-03  8.7629E-04  0.0000E+00 -2.6270E-03 -9.3530E-04
            -3.9131E-04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      862
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.6528E-04 -6.1542E-03 -2.6810E-04 -2.6306E-03 -1.9623E-02
 SE:             2.9816E-02  1.7267E-02  1.5957E-04  2.5158E-02  2.3921E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9290E-01  7.2154E-01  9.2926E-02  9.1672E-01  4.1203E-01

 ETASHRINKSD(%)  1.1285E-01  4.2152E+01  9.9465E+01  1.5717E+01  1.9863E+01
 ETASHRINKVR(%)  2.2558E-01  6.6536E+01  9.9997E+01  2.8964E+01  3.5780E+01
 EBVSHRINKSD(%)  4.5284E-01  4.2364E+01  9.9463E+01  1.5712E+01  1.8350E+01
 EBVSHRINKVR(%)  9.0363E-01  6.6781E+01  9.9997E+01  2.8956E+01  3.3333E+01
 RELATIVEINF(%)  9.8358E+01  8.3360E-01  2.9967E-04  2.2849E+00  5.0040E+00
 EPSSHRINKSD(%)  4.1804E+01
 EPSSHRINKVR(%)  6.6132E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1662.8497303976269     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -927.69890383388872     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.61
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1662.850       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  8.45E-01  9.92E-01  1.08E+00  8.95E-01  1.03E+00  1.04E+00  1.00E-02  8.49E-01  8.82E-01  1.08E+00
 


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
+        1.10E+03
 
 TH 2
+       -7.40E+00  4.94E+02
 
 TH 3
+        1.09E+01  1.76E+02  3.41E+02
 
 TH 4
+       -1.31E+01  4.99E+02 -1.28E+02  9.42E+02
 
 TH 5
+       -4.84E+00 -3.94E+02 -5.58E+02  1.49E+02  1.18E+03
 
 TH 6
+        3.24E-01 -1.44E+00  4.26E+00 -2.95E+00 -4.16E+00  1.83E+02
 
 TH 7
+        1.76E+00  1.59E+01  1.17E+01 -5.30E+00 -1.37E+01  1.90E-01  2.39E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.17E+00 -3.03E+01 -1.48E+01  2.23E+01  2.99E+00  2.11E+00  2.47E+01  0.00E+00  1.44E+02
 
 TH10
+       -3.98E-01  8.22E+00 -3.76E+01 -1.98E+01 -5.91E+01 -1.51E+00  1.34E+01  0.00E+00  9.45E-01  1.08E+02
 
 TH11
+       -6.70E+00 -1.89E+01 -3.06E+01 -7.49E+00  1.36E+01  2.21E-01  4.41E+00  0.00E+00  1.01E+01  2.89E+01  1.94E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.239
Stop Time:
Sat Sep 18 08:50:35 CDT 2021

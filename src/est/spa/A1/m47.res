Sat Sep 25 08:04:45 CDT 2021
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
$DATA ../../../../data/spa/A1/dat47.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1222.30627780133        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.4988E+01 -2.5083E+00  2.8925E+01 -5.9996E+01  1.1902E+02  3.3702E+01 -6.1983E+01 -1.3107E+01 -9.7277E+01 -1.1573E+02
            -5.9671E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1417.11962180265        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.7759E-01  8.2542E-01  9.1919E-01  1.1007E+00  8.6354E-01  8.7723E-01  1.2427E+00  9.6167E-01  1.3220E+00  1.2520E+00
             1.7949E+00
 PARAMETER:  7.7340E-02 -9.1868E-02  1.5739E-02  1.9596E-01 -4.6721E-02 -3.0983E-02  3.1726E-01  6.0914E-02  3.7918E-01  3.2475E-01
             6.8493E-01
 GRADIENT:  -8.7254E+00 -2.0794E+01 -1.8151E+01 -3.4211E+01  3.1248E+01 -1.1390E+01 -4.9226E+00  7.1992E+00  1.9124E+01  2.4277E+00
            -2.5308E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1427.92509827708        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.7200E-01  4.8259E-01  5.3611E-01  1.3934E+00  4.6858E-01  9.7574E-01  2.3708E+00  3.7465E-01  1.0859E+00  8.3676E-01
             1.6884E+00
 PARAMETER:  7.1603E-02 -6.2860E-01 -5.2342E-01  4.3175E-01 -6.5804E-01  7.5444E-02  9.6322E-01 -8.8176E-01  1.8244E-01 -7.8224E-02
             6.2380E-01
 GRADIENT:  -2.7162E+01  5.4881E+01  2.2092E+01  1.3346E+02 -3.0556E+01  2.2785E+01  1.7460E+01  1.3938E+00  3.2809E-02 -2.1897E-01
            -4.2432E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1445.55545746400        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.9179E-01  4.0574E-01  2.8235E-01  1.1823E+00  3.0696E-01  9.1587E-01  1.4539E+00  3.2631E-02  1.1026E+00  5.6992E-01
             1.8766E+00
 PARAMETER:  9.1760E-02 -8.0205E-01 -1.1646E+00  2.6742E-01 -1.0810E+00  1.2123E-02  4.7426E-01 -3.3225E+00  1.9763E-01 -4.6226E-01
             7.2948E-01
 GRADIENT:   3.6059E+00  1.0451E+01 -1.3452E+00  2.4273E+01  7.7177E-01 -1.4185E+00 -3.3802E+00 -3.4954E-04 -8.0267E+00 -2.6627E+00
             9.0426E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1445.86371973205        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.9014E-01  3.5456E-01  2.7512E-01  1.1790E+00  2.9203E-01  9.1936E-01  1.5660E+00  2.2873E-02  1.1228E+00  5.9299E-01
             1.8439E+00
 PARAMETER:  9.0089E-02 -9.3687E-01 -1.1906E+00  2.6463E-01 -1.1309E+00  1.5923E-02  5.4854E-01 -3.6778E+00  2.1586E-01 -4.2258E-01
             7.1186E-01
 GRADIENT:   1.3734E-01  3.7470E+00 -2.4092E+00  7.2910E+00  4.0685E+00 -1.5084E-01 -1.3637E+00  6.4589E-05 -2.9971E+00 -1.7899E+00
             2.8305E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1445.90323680231        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  9.9012E-01  3.1664E-01  2.7827E-01  1.1908E+00  2.8614E-01  9.1896E-01  1.6526E+00  1.8412E-02  1.1254E+00  6.2240E-01
             1.8383E+00
 PARAMETER:  9.0069E-02 -1.0500E+00 -1.1791E+00  2.7466E-01 -1.1513E+00  1.5489E-02  6.0236E-01 -3.8948E+00  2.1815E-01 -3.7417E-01
             7.0883E-01
 GRADIENT:  -5.9767E-03  2.7271E+00 -2.2339E-01  3.3212E+00  3.6536E-02 -2.9027E-02 -4.8989E-01  3.6074E-04 -1.3913E+00 -8.7415E-01
             7.2673E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1445.98866742936        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  9.8946E-01  2.5218E-01  2.8918E-01  1.2199E+00  2.8099E-01  9.1591E-01  1.8419E+00  1.1438E-02  1.1198E+00  6.7248E-01
             1.8468E+00
 PARAMETER:  8.9408E-02 -1.2776E+00 -1.1407E+00  2.9874E-01 -1.1694E+00  1.2159E-02  7.1078E-01 -4.3708E+00  2.1314E-01 -2.9678E-01
             7.1345E-01
 GRADIENT:   5.4849E-01  2.6121E+00  7.7700E+00  1.3407E+00 -1.4269E+01 -1.9177E-01 -1.2575E+00  4.9837E-04  3.2983E-01  2.3398E+00
             1.6220E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1448.72379784829        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      560
 NPARAMETR:  9.8459E-01  1.8937E-01  3.5605E-01  1.2833E+00  3.1805E-01  9.1578E-01  2.7330E+00  1.0000E-02  1.0997E+00  7.1507E-01
             1.8326E+00
 PARAMETER:  8.4473E-02 -1.5640E+00 -9.3270E-01  3.4944E-01 -1.0456E+00  1.2016E-02  1.1054E+00 -4.9047E+00  1.9500E-01 -2.3537E-01
             7.0575E-01
 GRADIENT:  -5.8311E+00  5.3491E+00  6.6063E+00 -3.4251E+01 -1.4776E+01  1.2543E+00  7.0927E+00  0.0000E+00  2.9461E+00  7.8143E-01
            -1.2263E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1449.70465153706        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      735
 NPARAMETR:  9.8635E-01  1.9550E-01  3.9602E-01  1.3366E+00  3.4646E-01  9.1214E-01  2.7049E+00  1.0000E-02  1.0637E+00  7.2470E-01
             1.8452E+00
 PARAMETER:  8.6257E-02 -1.5322E+00 -8.2628E-01  3.9014E-01 -9.6000E-01  8.0419E-03  1.0950E+00 -4.5329E+00  1.6179E-01 -2.2200E-01
             7.1257E-01
 GRADIENT:  -7.0128E-02  4.4603E-02 -2.6065E-02 -1.5722E-01  1.2628E-01  1.2494E-01 -6.6067E-03  0.0000E+00  6.4872E-02 -2.0282E-01
            -4.5557E-01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1449.70521249207        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      827
 NPARAMETR:  9.8634E-01  1.9435E-01  3.9640E-01  1.3374E+00  3.4656E-01  9.1177E-01  2.7200E+00  1.0000E-02  1.0630E+00  7.2545E-01
             1.8470E+00
 PARAMETER:  8.6243E-02 -1.5381E+00 -8.2534E-01  3.9074E-01 -9.5971E-01  7.6324E-03  1.1006E+00 -4.5528E+00  1.6112E-01 -2.2096E-01
             7.1355E-01
 GRADIENT:   5.5656E-03  1.1575E-02 -1.0069E-02 -2.4483E-02  9.8524E-03 -3.7145E-03  9.7550E-03  0.0000E+00  1.3554E-03  1.2813E-02
             1.2340E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      827
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.9214E-04  2.4590E-02 -1.9648E-04 -1.1434E-02  5.9383E-03
 SE:             2.9437E-02  1.2901E-02  2.4905E-04  2.7859E-02  2.2924E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8395E-01  5.6647E-02  4.3016E-01  6.8148E-01  7.9560E-01

 ETASHRINKSD(%)  1.3839E+00  5.6779E+01  9.9166E+01  6.6701E+00  2.3200E+01
 ETASHRINKVR(%)  2.7487E+00  8.1319E+01  9.9993E+01  1.2895E+01  4.1018E+01
 EBVSHRINKSD(%)  1.6036E+00  6.7169E+01  9.9144E+01  5.6073E+00  1.9193E+01
 EBVSHRINKVR(%)  3.1815E+00  8.9221E+01  9.9993E+01  1.0900E+01  3.4703E+01
 RELATIVEINF(%)  9.4941E+01  3.4768E+00  3.1772E-04  3.8875E+01  2.5497E+00
 EPSSHRINKSD(%)  3.9905E+01
 EPSSHRINKVR(%)  6.3886E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1449.7052124920651     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -714.55438592832695     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.51
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1449.705       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.94E-01  3.96E-01  1.34E+00  3.47E-01  9.12E-01  2.72E+00  1.00E-02  1.06E+00  7.25E-01  1.85E+00
 


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
+        1.33E+03
 
 TH 2
+       -5.83E+01  1.41E+03
 
 TH 3
+        7.73E+00  1.48E+02  5.68E+03
 
 TH 4
+       -1.63E+01  1.19E+02 -2.94E+02  4.97E+02
 
 TH 5
+        6.60E+01 -7.60E+02 -7.89E+03 -4.88E+01  1.23E+04
 
 TH 6
+        2.37E+00 -5.57E+00  9.91E+00 -4.75E+00  3.70E+00  2.24E+02
 
 TH 7
+        1.58E+00  9.83E+01 -4.85E+01 -5.50E+00  4.42E+01  1.86E-02  8.93E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.52E+00 -7.43E+01  8.69E+01  4.73E+00 -4.23E+01  2.04E+00 -4.90E+00  0.00E+00  1.39E+02
 
 TH10
+       -4.74E+00 -1.06E+02 -6.87E+01  5.35E+00 -5.25E+01  1.34E+00 -8.24E+00  0.00E+00  7.36E+00  1.83E+02
 
 TH11
+       -1.27E+01 -6.38E+00 -3.26E+01 -9.15E+00  3.56E+01  4.06E+00  5.04E-01  0.00E+00  1.16E+01  2.32E+01  6.58E+01
 
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
 #CPUT: Total CPU Time in Seconds,       15.367
Stop Time:
Sat Sep 25 08:05:02 CDT 2021

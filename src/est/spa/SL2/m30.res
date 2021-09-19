Sat Sep 18 12:11:24 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat30.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1672.77480680685        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.2551E+01 -1.7808E+01 -3.2078E+01  3.5545E+01  3.7511E+01  1.6687E+00 -1.6416E+01 -8.3785E-01  5.5951E+00 -1.3744E+01
            -1.1868E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1679.48286205439        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0145E+00  1.1822E+00  1.1152E+00  8.8165E-01  1.1415E+00  9.9718E-01  1.2537E+00  1.0414E+00  9.0995E-01  1.1612E+00
             1.0392E+00
 PARAMETER:  1.1439E-01  2.6735E-01  2.0907E-01 -2.5962E-02  2.3230E-01  9.7173E-02  3.2607E-01  1.4060E-01  5.6388E-03  2.4947E-01
             1.3844E-01
 GRADIENT:   5.7951E+01  1.3042E+01 -9.9756E-01  1.1157E+01  1.6832E+01  4.9843E-01  9.0958E+00 -5.7411E+00  3.8173E+00 -2.0307E+00
             4.4713E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1680.29067167783        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0093E+00  1.1239E+00  1.2682E+00  9.2752E-01  1.1693E+00  1.0189E+00  1.3027E+00  1.5052E+00  8.4743E-01  1.1658E+00
             1.0255E+00
 PARAMETER:  1.0926E-01  2.1682E-01  3.3759E-01  2.4763E-02  2.5640E-01  1.1873E-01  3.6446E-01  5.0891E-01 -6.5542E-02  2.5345E-01
             1.2513E-01
 GRADIENT:   4.7271E+01  1.0261E+01 -1.2127E+01  2.3648E+01  2.0875E+01  1.0693E+01  7.4773E+00  2.4386E+00  1.5061E+00 -1.8116E+00
             1.2893E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1680.98986967014        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      274
 NPARAMETR:  1.0145E+00  1.2194E+00  1.2400E+00  8.6489E-01  1.1802E+00  1.0156E+00  1.1603E+00  1.4803E+00  9.0945E-01  1.1842E+00
             1.0239E+00
 PARAMETER:  1.1440E-01  2.9834E-01  3.1512E-01 -4.5149E-02  2.6570E-01  1.1547E-01  2.4867E-01  4.9227E-01  5.0827E-03  2.6907E-01
             1.2361E-01
 GRADIENT:   5.4201E+00  5.8565E-01 -1.1299E+00  5.5918E+00  1.1771E+00  1.5182E+00 -1.4484E+00 -9.8127E-01 -4.3964E-01 -1.1661E-01
            -7.0568E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1681.31463557592        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      449
 NPARAMETR:  1.0127E+00  1.3863E+00  1.2110E+00  7.5902E-01  1.2387E+00  1.0150E+00  1.0706E+00  1.7567E+00  9.6796E-01  1.2130E+00
             1.0250E+00
 PARAMETER:  1.1261E-01  4.2664E-01  2.9147E-01 -1.7572E-01  3.1407E-01  1.1485E-01  1.6820E-01  6.6344E-01  6.7440E-02  2.9306E-01
             1.2472E-01
 GRADIENT:  -9.0714E-02  5.2117E+00  2.0548E+00  2.7128E+00 -2.9349E+00  8.3090E-01  1.6521E-01 -5.0448E-01 -3.0984E-01 -4.7541E-02
            -1.9595E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1681.32265166879        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      642             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0127E+00  1.3901E+00  1.2043E+00  7.5197E-01  1.2388E+00  1.0111E+00  1.0664E+00  1.7554E+00  9.7004E-01  1.2127E+00
             1.0259E+00
 PARAMETER:  1.1259E-01  4.2937E-01  2.8586E-01 -1.8506E-01  3.1412E-01  1.1104E-01  1.6429E-01  6.6268E-01  6.9578E-02  2.9283E-01
             1.2553E-01
 GRADIENT:   7.2152E+01 -4.5391E+01  1.1039E+02 -9.7224E+01 -1.8741E+01  3.0678E+00 -6.6942E+00 -6.2239E+01 -2.6274E+01 -1.9420E+00
            -1.9536E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1681.32695532168        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      799
 NPARAMETR:  1.0127E+00  1.3901E+00  1.2043E+00  7.5404E-01  1.2388E+00  1.0127E+00  1.0666E+00  1.7554E+00  9.7004E-01  1.2127E+00
             1.0255E+00
 PARAMETER:  1.1259E-01  4.2937E-01  2.8586E-01 -1.8231E-01  3.1412E-01  1.1265E-01  1.6448E-01  6.6268E-01  6.9578E-02  2.9283E-01
             1.2513E-01
 GRADIENT:  -1.6374E-01  2.2569E+00  2.3338E+00 -9.8396E-02 -3.3975E+00 -3.6324E-02 -2.0516E-02 -5.6603E-01 -3.8829E-01 -4.3236E-02
            -2.1707E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1681.34868665839        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      977
 NPARAMETR:  1.0141E+00  1.3904E+00  1.1783E+00  7.5151E-01  1.2380E+00  1.0139E+00  1.0648E+00  1.7285E+00  9.7050E-01  1.2126E+00
             1.0244E+00
 PARAMETER:  1.1398E-01  4.2961E-01  2.6405E-01 -1.8568E-01  3.1350E-01  1.1376E-01  1.6282E-01  6.4726E-01  7.0059E-02  2.9275E-01
             1.2413E-01
 GRADIENT:   2.7157E+00 -1.5665E+00  5.5239E-01 -5.8149E-01  2.0195E-01  3.7496E-01 -3.0431E-01 -1.4570E-01 -3.1559E-01  9.1296E-02
            -2.3918E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1681.35634814327        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1152
 NPARAMETR:  1.0132E+00  1.3906E+00  1.1505E+00  7.5109E-01  1.2296E+00  1.0131E+00  1.0657E+00  1.6840E+00  9.7458E-01  1.2051E+00
             1.0246E+00
 PARAMETER:  1.1308E-01  4.2976E-01  2.4017E-01 -1.8623E-01  3.0668E-01  1.1306E-01  1.6361E-01  6.2117E-01  7.4253E-02  2.8656E-01
             1.2433E-01
 GRADIENT:   6.5424E-01 -1.3287E+00 -4.0975E-02  1.9386E-01 -6.0435E-02  8.9934E-02  1.4719E-02  2.1512E-02  5.5371E-03 -5.0315E-02
             2.2053E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1681.35716342959        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1309
 NPARAMETR:  1.0128E+00  1.3913E+00  1.1521E+00  7.5090E-01  1.2304E+00  1.0129E+00  1.0655E+00  1.6870E+00  9.7471E-01  1.2061E+00
             1.0246E+00
 PARAMETER:  1.1270E-01  4.3026E-01  2.4159E-01 -1.8648E-01  3.0735E-01  1.1277E-01  1.6343E-01  6.2296E-01  7.4385E-02  2.8739E-01
             1.2435E-01
 GRADIENT:   5.3902E+01  3.1888E+01  1.4281E-01  1.0570E+01  1.7357E+00  7.3296E+00  1.3360E+00  1.3603E-01  4.3528E-01  4.2438E-01
             1.0374E-01

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1681.35727562726        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1437
 NPARAMETR:  1.0128E+00  1.3913E+00  1.1519E+00  7.5051E-01  1.2306E+00  1.0129E+00  1.0652E+00  1.6875E+00  9.7482E-01  1.2061E+00
             1.0246E+00
 PARAMETER:  1.1276E-01  4.3026E-01  2.4137E-01 -1.8701E-01  3.0747E-01  1.1283E-01  1.6319E-01  6.2323E-01  7.4494E-02  2.8743E-01
             1.2435E-01
 GRADIENT:  -1.7395E-02 -1.5452E+00 -1.3353E-02 -7.6473E-03  2.8376E-02 -9.9961E-04  5.9758E-03  8.7903E-03 -8.2021E-04 -2.0710E-03
             9.8837E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1437
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6182E-04 -1.0761E-02 -4.0425E-02  5.8986E-03 -4.2354E-02
 SE:             2.9853E-02  2.3629E-02  1.3540E-02  1.9525E-02  2.1841E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9567E-01  6.4881E-01  2.8312E-03  7.6257E-01  5.2474E-02

 ETASHRINKSD(%)  1.0000E-10  2.0840E+01  5.4638E+01  3.4589E+01  2.6830E+01
 ETASHRINKVR(%)  1.0000E-10  3.7337E+01  7.9423E+01  5.7214E+01  4.6462E+01
 EBVSHRINKSD(%)  4.3940E-01  2.0567E+01  6.0633E+01  3.7302E+01  2.2232E+01
 EBVSHRINKVR(%)  8.7686E-01  3.6905E+01  8.4502E+01  6.0690E+01  3.9522E+01
 RELATIVEINF(%)  9.8558E+01  1.8087E+00  1.4683E+00  1.0417E+00  1.9785E+01
 EPSSHRINKSD(%)  4.5031E+01
 EPSSHRINKVR(%)  6.9784E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1681.3572756272556     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -946.20644906351743     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.98
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1681.357       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.39E+00  1.15E+00  7.51E-01  1.23E+00  1.01E+00  1.07E+00  1.69E+00  9.75E-01  1.21E+00  1.02E+00
 


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
+        1.04E+03
 
 TH 2
+       -2.07E+08  3.33E+02
 
 TH 3
+        3.03E+00  2.78E+01  1.83E+08
 
 TH 4
+       -1.02E+01  4.14E+02 -8.54E+01  8.13E+02
 
 TH 5
+       -1.80E+00 -6.76E+01 -8.26E+01  1.08E+02  3.16E+02
 
 TH 6
+        3.96E-01 -2.07E+08  3.20E-01 -1.81E+00  2.17E-01  1.90E+02
 
 TH 7
+        2.73E+00  1.93E+01  2.16E+00 -1.84E+01 -3.30E+00  6.70E-01  7.61E+01
 
 TH 8
+       -6.06E-01 -8.46E+00 -4.83E+07  1.45E+01  3.77E-01 -1.04E-01  8.67E-01  6.97E+00
 
 TH 9
+        3.72E+00 -9.90E+00 -7.46E+00  6.98E+00  4.29E+00 -3.63E+00  2.87E+01  3.51E+00  3.48E+01
 
 TH10
+        3.75E-02 -3.98E+00 -6.21E+00  2.10E+00 -4.82E+01 -3.01E-03 -1.13E+00  3.88E+07  4.40E+00  5.98E+01
 
 TH11
+       -6.99E+00 -9.70E+00 -1.10E+01 -2.71E-01 -5.99E+00  6.22E-01  5.56E+00  4.38E+00  8.31E+00  8.92E+00  1.99E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.744
Stop Time:
Sat Sep 18 12:11:52 CDT 2021

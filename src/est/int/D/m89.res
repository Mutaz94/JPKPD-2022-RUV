Sat Sep 18 07:41:06 CDT 2021
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
$DATA ../../../../data/int/D/dat89.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   50565.4095516740        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.5900E+02  7.3542E+02  6.4503E+01  7.3468E+02  1.4952E+02 -3.8836E+03 -1.9027E+03 -1.3437E+02 -2.5105E+03 -1.1170E+03
            -9.8542E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -564.189891150173        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.2500E+00  1.8753E+00  8.2064E-01  2.0086E+00  1.0091E+00  4.4195E+00  3.5087E+00  9.9998E-01  1.6748E+00  1.3895E+00
             1.3012E+01
 PARAMETER:  3.2312E-01  7.2875E-01 -9.7670E-02  7.9746E-01  1.0909E-01  1.5860E+00  1.3553E+00  9.9983E-02  6.1571E-01  4.2893E-01
             2.6659E+00
 GRADIENT:  -1.0712E+01  3.1858E+01 -3.7734E+01  2.2439E+02  1.5813E+01  1.0041E+02 -8.1532E+01  3.8153E+00 -5.4931E+01  2.6363E+01
            -8.1281E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -624.684750349641        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.9344E-01  2.8821E+00  1.4434E+01  3.5770E+00  2.9803E+00  2.9984E+00  4.6580E+00  8.1616E-01  7.2090E+00  1.4704E+00
             1.2746E+01
 PARAMETER:  9.3417E-02  1.1585E+00  2.7696E+00  1.3745E+00  1.1920E+00  1.1981E+00  1.6386E+00 -1.0314E-01  2.0753E+00  4.8553E-01
             2.6452E+00
 GRADIENT:  -3.3831E+01  6.9220E+00 -7.2087E+00  5.4251E+01  3.5650E+01  8.3170E+01  5.0234E+00  8.4765E-01  1.9414E+01  2.5865E+01
             8.2069E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -696.055620299510        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0568E+00  2.1062E+00  1.4101E+01  2.2544E+00  2.0379E+00  2.2150E+00  5.5410E+00  3.2992E+00  2.0664E+00  8.1554E-01
             1.3577E+01
 PARAMETER:  1.5527E-01  8.4488E-01  2.7462E+00  9.1288E-01  8.1193E-01  8.9527E-01  1.8122E+00  1.2937E+00  8.2582E-01 -1.0391E-01
             2.7084E+00
 GRADIENT:  -4.9025E+01  4.4435E+01 -4.3030E+00  1.1931E+02 -5.2050E+01  2.0746E+01 -4.4150E+01  3.2439E+00 -2.3303E+01  5.9411E+00
             9.6235E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -749.311079605255        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0582E+00  1.2539E+00  1.7816E+01  1.3736E+00  2.4772E+00  2.0250E+00  5.3232E+00  4.9383E-01  1.5018E+00  5.0974E-01
             1.1720E+01
 PARAMETER:  1.5660E-01  3.2626E-01  2.9801E+00  4.1741E-01  1.0071E+00  8.0559E-01  1.7721E+00 -6.0557E-01  5.0667E-01 -5.7386E-01
             2.5613E+00
 GRADIENT:  -2.2435E+01  8.7421E-01 -3.2339E+00  2.3340E+00  1.2827E+01 -4.8156E+00 -1.5039E+01  1.8926E-02  1.5632E+00  2.8702E+00
            -3.5367E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -758.951335996838        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.1384E+00  4.9285E-01  7.9253E+01  2.0105E+00  2.5350E+00  2.0159E+00  7.6758E+00  1.0000E-02  1.7433E+00  3.6140E-02
             1.2159E+01
 PARAMETER:  2.2958E-01 -6.0756E-01  4.4726E+00  7.9840E-01  1.0302E+00  8.0107E-01  2.1381E+00 -5.0617E+00  6.5580E-01 -3.2204E+00
             2.5981E+00
 GRADIENT:   2.9996E+00  2.3980E+00 -7.5559E-01  2.9385E+01 -3.1804E+00 -2.1330E+00  5.7517E-01  0.0000E+00 -8.2218E-01  1.2003E-02
            -3.1868E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -759.218195911580        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  1.1297E+00  5.4711E-01  6.6676E+01  1.8702E+00  2.5256E+00  2.0212E+00  7.3241E+00  1.0136E-02  1.6612E+00  4.7742E-02
             1.2114E+01
 PARAMETER:  2.2196E-01 -5.0310E-01  4.2998E+00  7.2604E-01  1.0265E+00  8.0371E-01  2.0912E+00 -4.4917E+00  6.0754E-01 -2.9419E+00
             2.5944E+00
 GRADIENT:   8.9413E-01  5.6197E-01 -8.8135E-01  9.9602E+00  2.2554E-01 -1.3501E+00  1.2163E+00  1.4759E-06 -1.3963E+00  2.1072E-02
             1.8518E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -761.196654548075        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  1.1280E+00  2.2260E-01  2.7163E+02  2.0413E+00  2.5872E+00  2.0294E+00  9.7048E+00  1.0000E-02  1.7342E+00  1.0000E-02
             1.2088E+01
 PARAMETER:  2.2042E-01 -1.4024E+00  5.7044E+00  8.1361E-01  1.0506E+00  8.0775E-01  2.3726E+00 -7.9062E+00  6.5057E-01 -5.4744E+00
             2.5922E+00
 GRADIENT:   1.2965E+00  5.7242E-02 -2.5598E-01 -6.4391E+00  1.0961E+00  5.6382E-01  2.4800E+00  0.0000E+00 -3.2704E-01  0.0000E+00
             3.2969E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -761.443679997583        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  1.1175E+00  1.4963E-01  5.0135E+02  2.0790E+00  2.5982E+00  2.0166E+00  1.0514E+01  1.0000E-02  1.7081E+00  1.0000E-02
             1.2046E+01
 PARAMETER:  2.1111E-01 -1.7996E+00  6.3173E+00  8.3190E-01  1.0548E+00  8.0143E-01  2.4527E+00 -9.3390E+00  6.3538E-01 -6.4538E+00
             2.5887E+00
 GRADIENT:  -2.3557E+00 -6.5979E-01 -1.3335E-01  5.3480E+00 -1.0350E+00 -1.2249E+00 -2.7696E+00  0.0000E+00  1.0608E-01  0.0000E+00
            -6.3150E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -761.499697038014        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      975
 NPARAMETR:  1.1248E+00  1.3663E-01  5.9463E+02  2.0770E+00  2.6102E+00  2.0247E+00  1.0808E+01  1.0000E-02  1.7067E+00  1.0000E-02
             1.2088E+01
 PARAMETER:  2.1763E-01 -1.8905E+00  6.4879E+00  8.3094E-01  1.0594E+00  8.0542E-01  2.4803E+00 -9.6702E+00  6.3454E-01 -6.7662E+00
             2.5922E+00
 GRADIENT:  -1.7557E-01 -1.4542E-01 -1.4051E-01 -4.3638E-01  3.0307E-01 -1.5271E-02  2.9747E-01  0.0000E+00 -4.7932E-02  0.0000E+00
             5.1519E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -761.616777964082        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1068
 NPARAMETR:  1.1275E+00  1.4017E-01  3.9318E+04  2.0714E+00  2.6107E+00  2.0223E+00  1.0705E+01  1.0000E-02  1.6997E+00  1.0000E-02
             1.2053E+01
 PARAMETER:  2.2004E-01 -1.8649E+00  1.0679E+01  8.2823E-01  1.0596E+00  8.0422E-01  2.4707E+00 -9.6702E+00  6.3044E-01 -6.7662E+00
             2.5893E+00
 GRADIENT:  -1.5525E+01  1.2754E+00 -1.8802E-03  8.9476E+00 -9.5245E-01  8.9757E+00  2.7475E+01  0.0000E+00 -6.0216E+00  0.0000E+00
             6.1731E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -761.623725018232        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:     1269
 NPARAMETR:  1.1261E+00  1.4010E-01  3.9155E+04  2.0712E+00  2.6195E+00  2.0242E+00  1.0705E+01  1.0000E-02  1.7020E+00  1.0000E-02
             1.2064E+01
 PARAMETER:  2.1874E-01 -1.8654E+00  1.0675E+01  8.2814E-01  1.0630E+00  8.0516E-01  2.4707E+00 -9.6702E+00  6.3180E-01 -6.7662E+00
             2.5903E+00
 GRADIENT:  -8.1707E+00  8.5639E-01 -2.0081E-03  1.1536E+00  7.8237E-02  4.8652E+00  1.6338E+00  0.0000E+00 -4.9085E+00  0.0000E+00
             2.3370E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -761.627944560251        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1431
 NPARAMETR:  1.1239E+00  1.3980E-01  3.9421E+04  2.0717E+00  2.6187E+00  2.0235E+00  1.0716E+01  1.0000E-02  1.7054E+00  1.0000E-02
             1.2074E+01
 PARAMETER:  2.1683E-01 -1.8675E+00  1.0682E+01  8.2839E-01  1.0627E+00  8.0481E-01  2.4717E+00 -9.6702E+00  6.3380E-01 -6.7662E+00
             2.5911E+00
 GRADIENT:  -4.3052E+00  1.2901E-01 -1.9149E-03  5.6171E-01 -6.3293E-02  1.7488E+00 -8.1629E-02  0.0000E+00 -5.0492E-01  0.0000E+00
             1.5075E+00

0ITERATION NO.:   61    OBJECTIVE VALUE:  -761.627944560251        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     1461
 NPARAMETR:  1.1240E+00  1.3954E-01  3.9670E+04  2.0721E+00  2.6188E+00  2.0235E+00  1.0719E+01  1.0000E-02  1.7064E+00  1.0000E-02
             1.2068E+01
 PARAMETER:  2.1683E-01 -1.8675E+00  1.0682E+01  8.2839E-01  1.0627E+00  8.0481E-01  2.4717E+00 -9.6702E+00  6.3380E-01 -6.7662E+00
             2.5911E+00
 GRADIENT:  -2.6334E-02  2.0217E-01 -2.0241E-03 -2.5389E-01 -8.7956E-03 -8.1706E-03 -1.8160E-01  0.0000E+00 -2.2911E-01  0.0000E+00
             6.8313E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1461
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4314E-02  8.1149E-02  8.5787E-09 -7.9416E-02 -3.6490E-05
 SE:             2.8174E-02  1.9189E-02  5.0663E-09  1.9096E-02  9.6010E-05
 N:                     100         100         100         100         100

 P VAL.:         6.1142E-01  2.3502E-05  9.0402E-02  3.2010E-05  7.0390E-01

 ETASHRINKSD(%)  5.6132E+00  3.5713E+01  1.0000E+02  3.6027E+01  9.9678E+01
 ETASHRINKVR(%)  1.0911E+01  5.8672E+01  1.0000E+02  5.9074E+01  9.9999E+01
 EBVSHRINKSD(%)  7.5445E+00  4.5099E+01  1.0000E+02  2.2327E+01  9.9608E+01
 EBVSHRINKVR(%)  1.4520E+01  6.9858E+01  1.0000E+02  3.9669E+01  9.9998E+01
 RELATIVEINF(%)  8.5215E+01  2.3110E+01  0.0000E+00  4.4581E+01  3.4863E-04
 EPSSHRINKSD(%)  4.9308E+00
 EPSSHRINKVR(%)  9.6185E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -761.62794456025051     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       892.46141520816025     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    46.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    18.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -761.628       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.12E+00  1.40E-01  3.94E+04  2.07E+00  2.62E+00  2.02E+00  1.07E+01  1.00E-02  1.71E+00  1.00E-02  1.21E+01
 


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
+        6.26E+02
 
 TH 2
+       -1.46E+02  6.26E+03
 
 TH 3
+       -1.14E-04  9.66E-05  4.81E-11
 
 TH 4
+       -1.00E+02 -1.33E+03  2.99E-05  4.20E+02
 
 TH 5
+        1.24E+02  1.25E+02  7.61E-06 -2.65E+01  2.90E+01
 
 TH 6
+        1.90E+01  1.92E+02 -1.14E-05 -3.79E+01  5.62E+00  4.43E+01
 
 TH 7
+        1.48E+01  2.27E+02 -2.11E-06 -5.56E+01  3.53E+00  6.76E+00  9.61E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.51E+01 -7.10E+02  1.61E-05  1.94E+02 -1.45E+01 -3.49E+01 -2.77E+01  0.00E+00  1.65E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        8.29E-01  2.12E+02  3.30E-07 -5.47E+01  3.57E+00  9.48E+00  9.43E+00  0.00E+00 -2.57E+01  0.00E+00  1.97E+01
 
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
 #CPUT: Total CPU Time in Seconds,       65.356
Stop Time:
Sat Sep 18 07:42:13 CDT 2021

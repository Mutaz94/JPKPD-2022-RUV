Sat Sep 25 00:55:30 CDT 2021
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
$DATA ../../../../data/int/SL2/dat13.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      996
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      896
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
 RAW OUTPUT FILE (FILE): m13.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2476.87885011778        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3586E+01  4.8037E+01 -7.0990E+00 -4.4204E+01  5.1035E+01 -9.9924E-01 -8.6250E+01 -1.1848E+02 -3.2730E+01  3.2449E+00
            -2.5339E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3142.18163648858        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9610E-01  1.0701E+00  1.2395E+00  1.0054E+00  1.0976E+00  9.9993E-01  1.2618E+00  1.0714E+00  9.4544E-01  8.8447E-01
             1.8248E+00
 PARAMETER:  9.6088E-02  1.6779E-01  3.1475E-01  1.0539E-01  1.9311E-01  9.9933E-02  3.3251E-01  1.6900E-01  4.3893E-02 -2.2770E-02
             7.0148E-01
 GRADIENT:   3.8728E-01  7.0490E+00  7.7423E+00  1.6434E+00  1.0452E+01 -9.4192E-01  2.8874E+00 -9.7671E+00 -4.0750E+00  7.6357E-02
            -7.8884E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3143.29415061658        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0018E+00  1.0016E+00  1.2363E+00  1.0400E+00  1.0507E+00  1.0055E+00  1.3414E+00  1.3754E+00  9.1460E-01  8.1649E-01
             1.8064E+00
 PARAMETER:  1.0181E-01  1.0161E-01  3.1209E-01  1.3923E-01  1.4944E-01  1.0550E-01  3.9373E-01  4.1872E-01  1.0730E-02 -1.0275E-01
             6.9136E-01
 GRADIENT:   1.3935E+01 -7.4118E-01 -5.8511E-01  8.9576E-01  8.1617E+00  1.3220E+00  3.9316E+00 -2.6616E+00 -7.8484E+00  1.5048E+00
            -7.8408E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3144.79396904980        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9638E-01  9.9851E-01  1.1934E+00  1.0404E+00  1.0300E+00  1.0015E+00  1.3226E+00  1.3060E+00  9.4745E-01  7.8332E-01
             1.8709E+00
 PARAMETER:  9.6377E-02  9.8505E-02  2.7684E-01  1.3958E-01  1.2959E-01  1.0149E-01  3.7957E-01  3.6695E-01  4.6018E-02 -1.4421E-01
             7.2643E-01
 GRADIENT:  -4.3389E-01  8.2929E-01 -1.3621E-01  1.5320E-01  5.0004E-01 -1.8159E-01  8.0908E-02 -2.4321E-01 -2.6063E-01 -2.7376E-01
            -8.9961E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3144.79546461499        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      325
 NPARAMETR:  9.9651E-01  9.9559E-01  1.1927E+00  1.0416E+00  1.0280E+00  1.0018E+00  1.3250E+00  1.3093E+00  9.4765E-01  7.8166E-01
             1.8710E+00
 PARAMETER:  9.6500E-02  9.5583E-02  2.7624E-01  1.4077E-01  1.2759E-01  1.0177E-01  3.8138E-01  3.6946E-01  4.6235E-02 -1.4633E-01
             7.2649E-01
 GRADIENT:  -1.2621E+01 -1.7501E+00 -6.6253E-01 -4.1146E+00 -1.5854E+00 -1.3699E+00 -1.1598E+00 -1.9854E-01 -3.3946E-01 -2.0044E-01
            -1.5950E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3145.62880159332        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  1.0019E+00  1.0897E+00  1.7017E+00  1.0180E+00  1.2236E+00  1.0039E+00  1.2102E+00  2.2552E+00  8.6569E-01  9.6769E-01
             1.8543E+00
 PARAMETER:  1.0188E-01  1.8592E-01  6.3166E-01  1.1787E-01  3.0184E-01  1.0390E-01  2.9078E-01  9.1325E-01 -4.4225E-02  6.7151E-02
             7.1750E-01
 GRADIENT:  -9.7600E-01  5.9785E-01 -1.5555E+00 -1.6551E-01  2.1890E+00 -5.3381E-01 -8.7570E-01 -1.4300E-01  7.8664E-01  4.2347E+00
            -4.7194E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3146.06503948925        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  1.0019E+00  1.0242E+00  1.7413E+00  1.0590E+00  1.1811E+00  1.0044E+00  1.2945E+00  2.2909E+00  8.4149E-01  8.8764E-01
             1.8468E+00
 PARAMETER:  1.0187E-01  1.2389E-01  6.5464E-01  1.5736E-01  2.6643E-01  1.0442E-01  3.5815E-01  9.2893E-01 -7.2584E-02 -1.9191E-02
             7.1346E-01
 GRADIENT:  -5.3929E-01  2.3597E+00 -7.8111E-01  3.2030E+00 -2.1595E+00 -3.3555E-01  4.5432E-01 -1.3196E+00 -3.0282E-01  3.6787E+00
            -3.1788E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3146.07487344986        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      867
 NPARAMETR:  1.0021E+00  1.0242E+00  1.7410E+00  1.0591E+00  1.1810E+00  1.0053E+00  1.2946E+00  2.3248E+00  8.4188E-01  8.8762E-01
             1.8471E+00
 PARAMETER:  1.0212E-01  1.2392E-01  6.5448E-01  1.5740E-01  2.6637E-01  1.0528E-01  3.5824E-01  9.4362E-01 -7.2117E-02 -1.9216E-02
             7.1364E-01
 GRADIENT:  -1.2195E-02  2.4971E+00 -1.6156E+00  2.7102E+00 -2.6811E+00 -5.9065E-03  5.3006E-01  1.3822E-02 -7.8741E-03  3.6960E+00
            -1.2478E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3146.14880245137        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      959
 NPARAMETR:  9.9877E-01  1.0184E+00  1.7404E+00  1.0576E+00  1.1785E+00  1.0034E+00  1.2922E+00  2.3104E+00  8.4441E-01  8.5323E-01
             1.8484E+00
 PARAMETER:  9.8774E-02  1.1823E-01  6.5413E-01  1.5604E-01  2.6428E-01  1.0341E-01  3.5632E-01  9.3742E-01 -6.9115E-02 -5.8722E-02
             7.1433E-01
 GRADIENT:   5.6483E+00  1.3775E+00 -2.8400E-01  3.0545E-01  3.7441E+00  6.0256E-01 -8.8364E-02  3.4407E-01 -2.4519E-01  7.9774E-01
             2.5603E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3146.15099563671        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1029
 NPARAMETR:  9.9664E-01  1.0107E+00  1.7404E+00  1.0611E+00  1.1712E+00  1.0023E+00  1.3011E+00  2.3019E+00  8.4602E-01  8.3734E-01
             1.8483E+00
 PARAMETER:  9.6632E-02  1.1060E-01  6.5410E-01  1.5934E-01  2.5800E-01  1.0226E-01  3.6318E-01  9.3374E-01 -6.7210E-02 -7.7526E-02
             7.1425E-01
 GRADIENT:   1.0454E+00  1.2342E+00  4.1459E-01 -1.1350E+00  1.6254E+00  1.0688E-01 -1.9661E-01  4.2373E-02  1.9812E-03  3.9753E-01
             1.0900E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3146.16071442912        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1101
 NPARAMETR:  9.9478E-01  9.8887E-01  1.7403E+00  1.0727E+00  1.1546E+00  1.0013E+00  1.3276E+00  2.3033E+00  8.4399E-01  8.0284E-01
             1.8469E+00
 PARAMETER:  9.4762E-02  8.8809E-02  6.5406E-01  1.7018E-01  2.4375E-01  1.0134E-01  3.8334E-01  9.3434E-01 -6.9613E-02 -1.1961E-01
             7.1351E-01
 GRADIENT:  -2.8648E+00  5.2492E-01  6.3513E-01 -1.5385E+00 -7.8514E-01 -3.0441E-01 -1.6862E-01 -2.0499E-01  1.8443E-01 -8.1525E-02
             1.0268E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3146.25616527738        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1239
 NPARAMETR:  1.0022E+00  9.7584E-01  1.7403E+00  1.0839E+00  1.1494E+00  1.0053E+00  1.3541E+00  2.3394E+00  8.3763E-01  7.8776E-01
             1.8462E+00
 PARAMETER:  1.0216E-01  7.5542E-02  6.5409E-01  1.8055E-01  2.3925E-01  1.0533E-01  4.0313E-01  9.4990E-01 -7.7173E-02 -1.3856E-01
             7.1316E-01
 GRADIENT:   4.1521E-01 -3.2141E-01 -3.4206E+00  2.9708E-02  6.0232E-01  6.4450E-03  5.6337E-02 -2.8149E-02  2.3178E-02 -9.3335E-02
            -3.1399E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3146.28078652771        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1395
 NPARAMETR:  1.0016E+00  9.7774E-01  1.7674E+00  1.0827E+00  1.1504E+00  1.0051E+00  1.3505E+00  2.3388E+00  8.3785E-01  7.9097E-01
             1.8463E+00
 PARAMETER:  1.0164E-01  7.7487E-02  6.6950E-01  1.7945E-01  2.4012E-01  1.0512E-01  4.0051E-01  9.4962E-01 -7.6911E-02 -1.3449E-01
             7.1319E-01
 GRADIENT:   1.2411E+01  2.8741E+00  1.4743E+00  3.0180E+00 -1.9407E+00  1.3199E+00  1.0219E+00 -9.9525E-02  3.6113E-01  2.8993E-01
             9.9227E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3146.30978043348        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1553
 NPARAMETR:  1.0020E+00  9.7587E-01  1.7813E+00  1.0851E+00  1.1554E+00  1.0053E+00  1.3534E+00  2.3732E+00  8.3534E-01  7.8994E-01
             1.8461E+00
 PARAMETER:  1.0195E-01  7.5577E-02  6.7732E-01  1.8170E-01  2.4441E-01  1.0531E-01  4.0263E-01  9.6422E-01 -7.9918E-02 -1.3580E-01
             7.1308E-01
 GRADIENT:  -2.2174E-02  2.9024E-01 -1.0659E+00 -5.7486E-01 -1.5078E-01 -4.7780E-03 -5.8856E-02 -1.8953E-01  2.9695E-01 -6.4710E-02
             3.8571E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3146.31163721685        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1737
 NPARAMETR:  1.0020E+00  9.7472E-01  1.7819E+00  1.0861E+00  1.1554E+00  1.0053E+00  1.3548E+00  2.3785E+00  8.3293E-01  7.9071E-01
             1.8455E+00
 PARAMETER:  1.0197E-01  7.4392E-02  6.7769E-01  1.8262E-01  2.4442E-01  1.0533E-01  4.0365E-01  9.6645E-01 -8.2810E-02 -1.3483E-01
             7.1273E-01
 GRADIENT:   1.5988E-02  2.1409E-01 -1.3805E+00  2.7620E-02  3.7023E-01 -1.3531E-04 -2.5165E-02 -1.0586E-01 -4.4367E-03  2.5536E-02
            -1.3440E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -3146.31291244429        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1913
 NPARAMETR:  1.0019E+00  9.6995E-01  1.7819E+00  1.0887E+00  1.1521E+00  1.0053E+00  1.3607E+00  2.3832E+00  8.3198E-01  7.8426E-01
             1.8451E+00
 PARAMETER:  1.0195E-01  6.9493E-02  6.7769E-01  1.8502E-01  2.4156E-01  1.0532E-01  4.0797E-01  9.6843E-01 -8.3946E-02 -1.4302E-01
             7.1256E-01
 GRADIENT:   9.5010E-03  4.4925E-03 -1.5083E+00  4.2716E-04 -1.8424E-02  1.7106E-03 -1.2605E-03 -1.8880E-03 -5.0109E-04  1.2645E-03
             4.3167E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -3146.31591692968        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2075
 NPARAMETR:  1.0019E+00  9.6986E-01  1.7862E+00  1.0887E+00  1.1520E+00  1.0053E+00  1.3606E+00  2.3839E+00  8.3195E-01  7.8412E-01
             1.8451E+00
 PARAMETER:  1.0189E-01  6.9400E-02  6.8007E-01  1.8501E-01  2.4152E-01  1.0529E-01  4.0793E-01  9.6872E-01 -8.3986E-02 -1.4319E-01
             7.1256E-01
 GRADIENT:  -9.4767E-02  1.4343E-01 -1.0086E+00 -5.2058E-01 -8.1582E-01 -8.6296E-03 -3.7156E-02 -1.0290E-01  1.5604E-02  2.7265E-02
            -6.3641E-03

0ITERATION NO.:   85    OBJECTIVE VALUE:  -3146.31602332052        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2252
 NPARAMETR:  1.0019E+00  9.6962E-01  1.7859E+00  1.0889E+00  1.1523E+00  1.0053E+00  1.3609E+00  2.3866E+00  8.3182E-01  7.8373E-01
             1.8451E+00
 PARAMETER:  1.0194E-01  6.9153E-02  6.7990E-01  1.8520E-01  2.4176E-01  1.0531E-01  4.0813E-01  9.6987E-01 -8.4135E-02 -1.4369E-01
             7.1251E-01
 GRADIENT:  -1.8904E-03  1.9270E-02 -1.2490E+00 -2.6413E-01 -2.2307E-01 -2.3001E-03 -3.2434E-02 -6.2986E-03  2.7239E-02 -1.8334E-02
            -2.1983E-02

0ITERATION NO.:   88    OBJECTIVE VALUE:  -3146.31607295597        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     2344
 NPARAMETR:  1.0019E+00  9.6955E-01  1.7859E+00  1.0891E+00  1.1525E+00  1.0053E+00  1.3612E+00  2.3870E+00  8.3155E-01  7.8405E-01
             1.8451E+00
 PARAMETER:  1.0194E-01  6.9081E-02  6.7990E-01  1.8535E-01  2.4191E-01  1.0532E-01  4.0835E-01  9.7005E-01 -8.4460E-02 -1.4328E-01
             7.1252E-01
 GRADIENT:   3.6526E-03  9.9958E-03 -1.3373E+00 -1.5563E-02  2.1435E-02  1.3076E-03  3.3750E-03 -3.4920E-04  7.3722E-03  4.0014E-04
             2.0550E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2344
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2504E-03 -8.0633E-03 -2.7643E-02  4.0098E-03 -3.6467E-02
 SE:             2.9699E-02  2.3207E-02  2.2446E-02  2.3774E-02  1.8745E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6642E-01  7.2826E-01  2.1813E-01  8.6606E-01  5.1722E-02

 ETASHRINKSD(%)  5.0523E-01  2.2252E+01  2.4804E+01  2.0355E+01  3.7203E+01
 ETASHRINKVR(%)  1.0079E+00  3.9553E+01  4.3455E+01  3.6567E+01  6.0566E+01
 EBVSHRINKSD(%)  8.3140E-01  2.2397E+01  2.5686E+01  2.1731E+01  3.6870E+01
 EBVSHRINKVR(%)  1.6559E+00  3.9777E+01  4.4775E+01  3.8740E+01  6.0146E+01
 RELATIVEINF(%)  9.8331E+01  2.4946E+01  3.3543E+01  2.8437E+01  1.7591E+01
 EPSSHRINKSD(%)  1.9303E+01
 EPSSHRINKVR(%)  3.4879E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          896
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1646.7378515027735     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3146.3160729559722     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1499.5782214531987     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    62.30
 Elapsed covariance  time in seconds:    14.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3146.316       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  9.70E-01  1.79E+00  1.09E+00  1.15E+00  1.01E+00  1.36E+00  2.39E+00  8.32E-01  7.84E-01  1.85E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.08E-02  2.91E-01  2.76E-03  1.64E-01  1.88E-01  7.08E-02  4.03E-01  3.93E-01  1.07E-01  3.86E-01  2.12E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.46E-04
 
 TH 2
+       -7.57E-04  8.48E-02
 
 TH 3
+       -9.36E-06  4.19E-04  7.62E-06
 
 TH 4
+        4.66E-04 -4.67E-02 -2.46E-04  2.69E-02
 
 TH 5
+       -4.79E-04  5.29E-02  3.04E-04 -2.94E-02  3.55E-02
 
 TH 6
+       -1.23E-04  3.12E-03  1.21E-05 -1.37E-03  1.52E-03  5.01E-03
 
 TH 7
+        1.51E-03 -1.10E-01 -5.56E-04  6.05E-02 -7.04E-02 -3.30E-03  1.63E-01
 
 TH 8
+       -5.07E-04 -5.32E-02 -7.19E-04  3.12E-02 -3.68E-02  9.87E-04  6.60E-02  1.54E-01
 
 TH 9
+        2.05E-04  1.57E-02  1.33E-04 -8.84E-03  9.39E-03  1.27E-03 -1.97E-02 -1.76E-02  1.14E-02
 
 TH10
+       -2.30E-03  1.00E-01  4.25E-04 -5.50E-02  6.61E-02  4.23E-03 -1.38E-01 -4.83E-02  1.33E-02  1.49E-01
 
 TH11
+       -1.86E-04  3.79E-02  1.07E-04 -2.15E-02  2.21E-02 -1.34E-04 -4.81E-02 -3.47E-02  6.23E-03  3.49E-02  4.50E-02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.08E-02
 
 TH 2
+       -8.46E-02  2.91E-01
 
 TH 3
+       -1.10E-01  5.21E-01  2.76E-03
 
 TH 4
+        9.23E-02 -9.77E-01 -5.45E-01  1.64E-01
 
 TH 5
+       -8.26E-02  9.65E-01  5.84E-01 -9.52E-01  1.88E-01
 
 TH 6
+       -5.66E-02  1.52E-01  6.18E-02 -1.18E-01  1.14E-01  7.08E-02
 
 TH 7
+        1.21E-01 -9.37E-01 -4.99E-01  9.15E-01 -9.26E-01 -1.16E-01  4.03E-01
 
 TH 8
+       -4.20E-02 -4.66E-01 -6.64E-01  4.84E-01 -4.98E-01  3.55E-02  4.17E-01  3.93E-01
 
 TH 9
+        6.24E-02  5.06E-01  4.53E-01 -5.06E-01  4.68E-01  1.68E-01 -4.57E-01 -4.21E-01  1.07E-01
 
 TH10
+       -1.94E-01  8.92E-01  3.99E-01 -8.67E-01  9.08E-01  1.55E-01 -8.86E-01 -3.19E-01  3.23E-01  3.86E-01
 
 TH11
+       -2.85E-02  6.14E-01  1.83E-01 -6.19E-01  5.54E-01 -8.90E-03 -5.62E-01 -4.17E-01  2.76E-01  4.26E-01  2.12E-01
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.27E+03
 
 TH 2
+       -4.55E+01  4.87E+02
 
 TH 3
+        6.12E+03  1.71E+03  4.08E+05
 
 TH 4
+       -4.84E+01  4.42E+02  2.56E+03  9.19E+02
 
 TH 5
+       -1.95E+02 -2.18E+02 -6.29E+03  7.68E+01  7.02E+02
 
 TH 6
+        1.09E+01 -6.39E+01 -7.11E+02 -5.49E+01  4.70E+01  2.23E+02
 
 TH 7
+       -8.75E+00  4.90E+01  4.84E+02  1.35E+01  1.90E+01 -1.07E+01  5.91E+01
 
 TH 8
+        1.86E+01 -2.51E+00  1.30E+03  2.93E+00  6.99E+00 -5.73E+00  3.10E+00  1.49E+01
 
 TH 9
+       -2.44E+01 -4.57E+01 -5.45E+02  2.24E+01  2.08E+01 -2.49E+01  7.65E+00  6.70E+00  1.48E+02
 
 TH10
+        9.21E+01 -1.71E+01  1.98E+03  4.76E+00 -1.03E+02 -1.28E+01  1.59E+01  1.17E+00  2.44E+01  6.88E+01
 
 TH11
+        3.89E+01 -2.62E+01  1.98E+03  3.32E+01 -7.82E+00  4.38E+00  6.86E+00  1.00E+01  1.40E+01  2.38E+01  5.42E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       76.603
Stop Time:
Sat Sep 25 00:56:48 CDT 2021

Sat Sep 25 14:28:34 CDT 2021
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
$DATA ../../../../data/spa/D/dat59.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m59.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   11834.3130174687        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.9962E+02  2.5315E+02 -4.6553E+01  1.7755E+01  4.2937E+02 -2.9606E+03 -7.6441E+02 -4.7806E+01 -1.5468E+03 -1.0238E+03
            -2.0886E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -675.897318007054        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.7205E+00  6.3107E-01  7.7378E-01  2.0982E+00  2.5228E+00  3.0971E+00  1.7139E+00  9.2334E-01  2.0597E+00  1.2630E+00
             1.1953E+01
 PARAMETER:  6.4260E-01 -3.6034E-01 -1.5647E-01  8.4106E-01  1.0254E+00  1.2305E+00  6.3879E-01  2.0244E-02  8.2257E-01  3.3345E-01
             2.5810E+00
 GRADIENT:   3.5931E+01  3.2832E+01 -3.2574E+01  3.6284E+01 -4.4481E+01  5.3878E+01  1.9718E+00  3.1891E+00 -1.4779E+01 -2.0037E+00
             7.8998E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -698.653409085559        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.5832E+00  5.5335E-01  7.9694E-01  2.1665E+00  4.6796E+00  3.3297E+00  8.2928E+00  5.1917E-01  2.2304E+00  8.6317E-01
             1.0880E+01
 PARAMETER:  5.5943E-01 -4.9177E-01 -1.2698E-01  8.7312E-01  1.6432E+00  1.3029E+00  2.2154E+00 -5.5553E-01  9.0218E-01 -4.7139E-02
             2.4869E+00
 GRADIENT:   3.2304E+01  2.6910E+01 -1.4509E+01  5.3453E+01 -1.3223E+01  6.6443E+01  6.1571E+00 -1.1046E+00  7.9107E+00  3.7410E-01
             8.2897E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -779.614040082583        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.3327E-01  2.3629E-01  3.8769E-01  1.3361E+00  7.4216E+00  2.3690E+00  9.1051E+00  2.9857E+00  1.1993E+00  2.1520E+00
             1.0015E+01
 PARAMETER:  3.0935E-02 -1.3427E+00 -8.4754E-01  3.8978E-01  2.1044E+00  9.6246E-01  2.3088E+00  1.1938E+00  2.8175E-01  8.6639E-01
             2.4041E+00
 GRADIENT:  -5.8139E+01  3.8417E+01 -3.2506E+01  4.6112E+01 -1.3129E+01  3.0301E+01  6.1354E+00  1.2529E+01  9.3156E+00  3.7776E-01
             8.8795E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -860.158640608951        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  6.5971E-01  1.6597E-02  5.1994E-02  4.5895E-01  3.3754E+01  1.9321E+00  4.1069E+00  1.7544E+00  2.6326E-01  8.4111E+00
             8.1084E+00
 PARAMETER: -3.1595E-01 -3.9985E+00 -2.8566E+00 -6.7881E-01  3.6191E+00  7.5861E-01  1.5127E+00  6.6210E-01 -1.2346E+00  2.2296E+00
             2.1929E+00
 GRADIENT:   3.3447E+01 -6.6910E-01  5.8447E+00 -1.4882E+01  7.3435E-02 -5.2221E+00  2.0270E-02  1.8190E+01  1.4692E+00  3.0711E-04
            -1.9941E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -862.303715065947        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  5.6482E-01  1.0963E-02  3.8192E-02  3.7889E-01  4.9933E+01  1.9143E+00  2.4577E+00  1.4679E+00  2.1084E-01  7.5254E+00
             8.2271E+00
 PARAMETER: -4.7125E-01 -4.4132E+00 -3.1651E+00 -8.7050E-01  4.0107E+00  7.4937E-01  9.9922E-01  4.8382E-01 -1.4567E+00  2.1183E+00
             2.2074E+00
 GRADIENT:   2.3512E-01 -3.7667E-02 -1.1907E+01  3.0745E+01  3.6810E-02 -7.6973E+00  3.8721E-03 -1.0998E+00  1.9746E-01 -3.0623E-04
            -7.3415E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -862.305542878403        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      456
 NPARAMETR:  5.5966E-01  1.0611E-02  3.7288E-02  3.7223E-01  5.1347E+01  1.9144E+00  2.3963E+00  1.4591E+00  2.0663E-01  7.4975E+00
             8.2246E+00
 PARAMETER: -4.8043E-01 -4.4459E+00 -3.1891E+00 -8.8823E-01  4.0386E+00  7.4940E-01  9.7391E-01  4.7785E-01 -1.4768E+00  2.1146E+00
             2.2071E+00
 GRADIENT:   1.5491E-01 -3.4049E-02 -1.1706E+01  3.0228E+01  3.4360E-02 -7.5695E+00  3.5505E-03 -1.1590E+00  1.7844E-01 -2.4435E-04
            -7.1848E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -862.720647400724        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      618
 NPARAMETR:  5.7259E-01  1.1404E-02  3.9519E-02  3.8144E-01  4.9248E+01  1.9628E+00  2.2375E+00  1.4728E+00  2.1351E-01  7.0899E+00
             8.3145E+00
 PARAMETER: -4.5758E-01 -4.3738E+00 -3.1310E+00 -8.6380E-01  3.9969E+00  7.7436E-01  9.0535E-01  4.8714E-01 -1.4441E+00  2.0587E+00
             2.2180E+00
 GRADIENT:  -2.8000E-01 -9.2789E-02  1.1116E-02  1.7247E-02  3.0278E-02 -9.3585E-02  3.3548E-03  3.1061E-02  3.0671E-01 -1.7095E-04
             9.4423E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -862.751835533234        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  5.7440E-01  1.2032E-02  4.0055E-02  3.8466E-01  4.5316E+01  1.9659E+00  2.1194E+00  1.4860E+00  1.9086E-01  1.1725E+01
             8.2677E+00
 PARAMETER: -4.5444E-01 -4.3202E+00 -3.1175E+00 -8.5541E-01  3.9137E+00  7.7597E-01  8.5112E-01  4.9611E-01 -1.5562E+00  2.5617E+00
             2.2124E+00
 GRADIENT:  -5.2099E-01 -1.1461E-01  5.7070E-01 -3.8087E-01  3.3481E-02  2.1959E-01  3.1510E-03  6.2163E-02  2.3981E-01  1.5972E-04
            -3.0173E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -862.756098804511        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  5.7461E-01  1.2103E-02  4.0088E-02  3.8488E-01  4.4961E+01  1.9659E+00  2.0971E+00  1.4869E+00  1.8838E-01  1.2343E+01
             8.2686E+00
 PARAMETER: -4.5406E-01 -4.3143E+00 -3.1167E+00 -8.5482E-01  3.9058E+00  7.7596E-01  8.4057E-01  4.9670E-01 -1.5693E+00  2.6131E+00
             2.2125E+00
 GRADIENT:  -5.0996E-01 -1.1818E-01  5.6428E-01 -3.6539E-01  3.3744E-02  2.1811E-01  3.0954E-03  6.1148E-02  2.3280E-01  3.9525E-04
            -2.9567E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -862.758519791023        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1164
 NPARAMETR:  5.7473E-01  1.2146E-02  4.0106E-02  3.8500E-01  4.4758E+01  1.9659E+00  2.0841E+00  1.4874E+00  1.8701E-01  1.2712E+01
             8.2691E+00
 PARAMETER: -4.5385E-01 -4.3108E+00 -3.1162E+00 -8.5450E-01  3.9013E+00  7.7595E-01  8.3432E-01  4.9703E-01 -1.5766E+00  2.6425E+00
             2.2125E+00
 GRADIENT:  -5.0099E-01 -1.2048E-01  5.6050E-01 -3.5776E-01  3.3870E-02  2.1529E-01  3.0645E-03  6.0240E-02  2.2896E-01  5.8367E-04
            -2.9175E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -862.839662151444        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1333
 NPARAMETR:  5.7204E-01  1.2425E-02  4.0238E-02  3.8471E-01  4.4260E+01  1.9531E+00  2.0134E+00  1.5018E+00  5.5241E-02  1.2793E+01
             8.2882E+00
 PARAMETER: -4.5855E-01 -4.2880E+00 -3.1129E+00 -8.5528E-01  3.8901E+00  7.6941E-01  7.9984E-01  5.0667E-01 -2.7960E+00  2.6489E+00
             2.2148E+00
 GRADIENT:  -3.4837E+00 -1.2905E-01  2.8919E+00 -2.0064E+00  3.0430E-02 -2.1304E+00  2.7276E-03 -6.8701E-01  1.5386E-02  5.5624E-04
            -1.0532E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -862.891415484895        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1517
 NPARAMETR:  5.8205E-01  1.5212E-02  4.0233E-02  3.8471E-01  4.1235E+01  1.9804E+00  1.6346E+00  1.5052E+00  2.5564E-02  1.2792E+01
             8.2893E+00
 PARAMETER: -4.4119E-01 -4.0857E+00 -3.1131E+00 -8.5527E-01  3.8193E+00  7.8329E-01  5.9137E-01  5.0891E-01 -3.5666E+00  2.6488E+00
             2.2150E+00
 GRADIENT:   5.2866E+00 -5.9446E-01  2.3259E+00 -5.3098E+00  5.2634E-02  3.1064E+00  8.2956E-03 -4.2999E-01  3.5603E-03  1.1558E-03
            -1.8697E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -862.915860935247        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1618
 NPARAMETR:  5.7614E-01  1.5228E-02  4.0264E-02  3.8463E-01  4.1274E+01  1.9802E+00  1.6120E+00  1.5054E+00  2.1852E-02  1.2784E+01
             8.2939E+00
 PARAMETER: -4.5141E-01 -4.0846E+00 -3.1123E+00 -8.5548E-01  3.8202E+00  7.8320E-01  5.7749E-01  5.0904E-01 -3.7235E+00  2.6482E+00
             2.2155E+00
 GRADIENT:   4.1781E+00 -5.5696E-01  9.1764E+00 -8.5586E-01  4.6645E-02  5.5597E+00  8.4758E-03 -2.6255E-01  3.3592E-03  1.9967E-03
             1.3623E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -862.915868021700        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1693
 NPARAMETR:  5.7593E-01  1.5228E-02  4.0263E-02  3.8463E-01  4.1274E+01  1.9802E+00  1.6025E+00  1.5054E+00  2.1369E-02  1.2784E+01
             8.2939E+00
 PARAMETER: -4.5176E-01 -4.0846E+00 -3.1123E+00 -8.5548E-01  3.8202E+00  7.8319E-01  5.7156E-01  5.0904E-01 -3.7458E+00  2.6482E+00
             2.2155E+00
 GRADIENT:   3.9898E+00 -5.5612E-01  9.1803E+00 -7.6685E-01  4.6458E-02  5.5553E+00  8.3988E-03 -2.6756E-01  3.2144E-03  2.0003E-03
             1.3930E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -862.915873478043        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1768
 NPARAMETR:  5.7575E-01  1.5228E-02  4.0262E-02  3.8463E-01  4.1274E+01  1.9802E+00  1.5884E+00  1.5054E+00  2.0924E-02  1.2784E+01
             8.2939E+00
 PARAMETER: -4.5208E-01 -4.0846E+00 -3.1123E+00 -8.5548E-01  3.8202E+00  7.8319E-01  5.6274E-01  5.0904E-01 -3.7669E+00  2.6482E+00
             2.2155E+00
 GRADIENT:   3.8237E+00 -5.5582E-01  9.1790E+00 -6.8119E-01  4.6307E-02  5.5510E+00  8.2808E-03 -2.7212E-01  3.0838E-03  2.0038E-03
             1.4186E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -862.915878329347        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1843
 NPARAMETR:  5.7559E-01  1.5228E-02  4.0262E-02  3.8463E-01  4.1274E+01  1.9802E+00  1.5715E+00  1.5054E+00  2.0512E-02  1.2784E+01
             8.2939E+00
 PARAMETER: -4.5236E-01 -4.0846E+00 -3.1124E+00 -8.5548E-01  3.8202E+00  7.8319E-01  5.5206E-01  5.0904E-01 -3.7867E+00  2.6482E+00
             2.2155E+00
 GRADIENT:   3.6740E+00 -5.5588E-01  9.1734E+00 -5.9742E-01  4.6181E-02  5.5468E+00  8.1385E-03 -2.7628E-01  2.9644E-03  2.0070E-03
             1.4403E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -862.915882420335        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1918
 NPARAMETR:  5.7544E-01  1.5228E-02  4.0261E-02  3.8463E-01  4.1274E+01  1.9802E+00  1.5526E+00  1.5054E+00  2.0127E-02  1.2784E+01
             8.2938E+00
 PARAMETER: -4.5262E-01 -4.0846E+00 -3.1124E+00 -8.5548E-01  3.8202E+00  7.8319E-01  5.3993E-01  5.0904E-01 -3.8057E+00  2.6482E+00
             2.2155E+00
 GRADIENT:   3.5354E+00 -5.5619E-01  9.1642E+00 -5.1378E-01  4.6074E-02  5.5425E+00  7.9789E-03 -2.8021E-01  2.8559E-03  2.0101E-03
             1.4591E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -862.915886059444        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1993
 NPARAMETR:  5.7529E-01  1.5228E-02  4.0259E-02  3.8463E-01  4.1275E+01  1.9802E+00  1.5319E+00  1.5054E+00  1.9748E-02  1.2784E+01
             8.2938E+00
 PARAMETER: -4.5287E-01 -4.0846E+00 -3.1124E+00 -8.5549E-01  3.8202E+00  7.8319E-01  5.2651E-01  5.0904E-01 -3.8247E+00  2.6482E+00
             2.2155E+00
 GRADIENT:   3.4031E+00 -5.5669E-01  9.1511E+00 -4.2786E-01  4.5980E-02  5.5382E+00  7.8052E-03 -2.8396E-01  2.7514E-03  2.0131E-03
             1.4758E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -862.915889212232        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     2068
 NPARAMETR:  5.7516E-01  1.5228E-02  4.0258E-02  3.8462E-01  4.1275E+01  1.9802E+00  1.5106E+00  1.5054E+00  1.9392E-02  1.2784E+01
             8.2937E+00
 PARAMETER: -4.5311E-01 -4.0846E+00 -3.1125E+00 -8.5549E-01  3.8203E+00  7.8319E-01  5.1251E-01  5.0904E-01 -3.8429E+00  2.6482E+00
             2.2155E+00
 GRADIENT:   3.2804E+00 -5.5731E-01  9.1350E+00 -3.4244E-01  4.5899E-02  5.5339E+00  7.6278E-03 -2.8746E-01  2.6538E-03  2.0158E-03
             1.4901E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -862.915892006938        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     2143
 NPARAMETR:  5.7503E-01  1.5229E-02  4.0256E-02  3.8462E-01  4.1275E+01  1.9802E+00  1.4883E+00  1.5054E+00  1.9041E-02  1.2783E+01
             8.2937E+00
 PARAMETER: -4.5334E-01 -4.0846E+00 -3.1125E+00 -8.5549E-01  3.8203E+00  7.8319E-01  4.9761E-01  5.0904E-01 -3.8612E+00  2.6482E+00
             2.2155E+00
 GRADIENT:   3.1615E+00 -5.5805E-01  9.1156E+00 -2.5411E-01  4.5828E-02  5.5294E+00  7.4419E-03 -2.9088E-01  2.5593E-03  2.0185E-03
             1.5026E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -862.915893499971        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     2219
 NPARAMETR:  5.7491E-01  1.5229E-02  4.0254E-02  3.8462E-01  4.1275E+01  1.9802E+00  1.4668E+00  1.5054E+00  1.8717E-02  1.2783E+01
             8.2936E+00
 PARAMETER: -4.5354E-01 -4.0846E+00 -3.1125E+00 -8.5549E-01  3.8203E+00  7.8319E-01  4.8307E-01  5.0904E-01 -3.8783E+00  2.6481E+00
             2.2155E+00
 GRADIENT:   3.0539E+00 -5.5880E-01  9.0946E+00 -1.6920E-01  4.5768E-02  5.5252E+00  7.2647E-03 -2.9398E-01  2.4740E-03  2.0208E-03
             1.5129E+00

0ITERATION NO.:  110    OBJECTIVE VALUE:  -862.915896121601        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     2295
 NPARAMETR:  5.7478E-01  1.5229E-02  4.0252E-02  3.8462E-01  4.1275E+01  1.9802E+00  1.4428E+00  1.5054E+00  1.8368E-02  1.2783E+01
             8.2935E+00
 PARAMETER: -4.5376E-01 -4.0845E+00 -3.1126E+00 -8.5549E-01  3.8203E+00  7.8319E-01  4.6658E-01  5.0904E-01 -3.8971E+00  2.6481E+00
             2.2155E+00
 GRADIENT:   2.9404E+00 -5.5970E-01  9.0682E+00 -7.3595E-02  4.5711E-02  5.5205E+00  7.0678E-03 -2.9727E-01  2.3844E-03  2.0232E-03
             1.5225E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -862.915897799354        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:     2372
 NPARAMETR:  5.7468E-01  1.5230E-02  4.0250E-02  3.8462E-01  4.1276E+01  1.9802E+00  1.4214E+00  1.5054E+00  1.8067E-02  1.2783E+01
             8.2935E+00
 PARAMETER: -4.5395E-01 -4.0845E+00 -3.1126E+00 -8.5550E-01  3.8203E+00  7.8319E-01  4.5163E-01  5.0904E-01 -3.9137E+00  2.6481E+00
             2.2155E+00
 GRADIENT:   2.8432E+00 -5.6053E-01  9.0422E+00  1.3108E-02  4.5667E-02  5.5163E+00  6.8940E-03 -3.0006E-01  2.3070E-03  2.0251E-03
             1.5296E+00

0ITERATION NO.:  120    OBJECTIVE VALUE:  -862.915898536730        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:     2450
 NPARAMETR:  5.7457E-01  1.5230E-02  4.0248E-02  3.8462E-01  4.1276E+01  1.9802E+00  1.4007E+00  1.5054E+00  1.7784E-02  1.2783E+01
             8.2934E+00
 PARAMETER: -4.5413E-01 -4.0845E+00 -3.1127E+00 -8.5550E-01  3.8203E+00  7.8319E-01  4.3699E-01  5.0904E-01 -3.9295E+00  2.6481E+00
             2.2155E+00
 GRADIENT:   2.7524E+00 -5.6135E-01  9.0150E+00  9.8308E-02  4.5629E-02  5.5121E+00  6.7267E-03 -3.0267E-01  2.2365E-03  2.0269E-03
             1.5352E+00

0ITERATION NO.:  125    OBJECTIVE VALUE:  -862.915899377303        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     2526
 NPARAMETR:  5.7445E-01  1.5230E-02  4.0246E-02  3.8462E-01  4.1276E+01  1.9802E+00  1.3748E+00  1.5054E+00  1.7434E-02  1.2783E+01
             8.2933E+00
 PARAMETER: -4.5434E-01 -4.0845E+00 -3.1128E+00 -8.5550E-01  3.8203E+00  7.8318E-01  4.1828E-01  5.0904E-01 -3.9493E+00  2.6481E+00
             2.2154E+00
 GRADIENT:   2.6418E+00 -5.6240E-01  8.9777E+00  2.0799E-01  4.5588E-02  5.5068E+00  6.5186E-03 -3.0584E-01  2.1498E-03  2.0289E-03
             1.5409E+00

0ITERATION NO.:  130    OBJECTIVE VALUE:  -862.926816082053        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2712
 NPARAMETR:  5.7563E-01  1.5243E-02  4.0178E-02  3.8461E-01  4.1284E+01  1.9801E+00  5.5391E-01  1.5054E+00  1.0000E-02  1.2781E+01
             8.2920E+00
 PARAMETER: -4.5228E-01 -4.0836E+00 -3.1144E+00 -8.5554E-01  3.8205E+00  7.8314E-01 -4.9076E-01  5.0908E-01 -9.3617E+00  2.6480E+00
             2.2153E+00
 GRADIENT:  -6.5099E-01 -6.1653E-01  2.0194E+00 -1.9222E+00  4.8081E-02  2.9690E+00  1.2176E-03 -5.1693E-01  0.0000E+00  1.3203E-03
            -6.7532E-01

0ITERATION NO.:  132    OBJECTIVE VALUE:  -862.926867211843        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:     2789
 NPARAMETR:  5.7563E-01  1.5246E-02  4.0183E-02  3.8461E-01  4.1292E+01  1.9801E+00  5.4931E-01  1.5055E+00  1.0000E-02  1.2780E+01
             8.2928E+00
 PARAMETER: -4.5228E-01 -4.0836E+00 -3.1145E+00 -8.5554E-01  3.8205E+00  7.8314E-01 -4.9859E-01  5.0908E-01 -9.4041E+00  2.6480E+00
             2.2153E+00
 GRADIENT:   2.2720E+03 -1.2714E+02 -3.6099E+02 -2.0112E+03 -2.7703E+02  2.2892E+03  1.3304E-03 -1.0272E+03  0.0000E+00  2.0497E+02
            -4.8356E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2789
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.9451E-03 -1.5755E-04  6.3902E-03 -3.4108E-04 -9.4205E-04
 SE:             2.8987E-02  6.7506E-05  2.4865E-02  2.4619E-04  3.9236E-04
 N:                     100         100         100         100         100

 P VAL.:         8.3750E-01  1.9601E-02  7.9718E-01  1.6592E-01  1.6352E-02

 ETASHRINKSD(%)  2.8896E+00  9.9774E+01  1.6699E+01  9.9175E+01  9.8686E+01
 ETASHRINKVR(%)  5.6956E+00  9.9999E+01  3.0609E+01  9.9993E+01  9.9983E+01
 EBVSHRINKSD(%)  1.9359E+00  9.9666E+01  1.5975E+01  9.9181E+01  9.8898E+01
 EBVSHRINKVR(%)  3.8344E+00  9.9999E+01  2.9397E+01  9.9993E+01  9.9988E+01
 RELATIVEINF(%)  5.6901E+00  4.6132E-04  1.2444E+00  8.0844E-05  1.4161E-03
 EPSSHRINKSD(%)  1.5103E+01
 EPSSHRINKVR(%)  2.7924E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -862.92686721184339     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -127.77604064810521     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.30
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -862.927       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.76E-01  1.52E-02  4.02E-02  3.85E-01  4.13E+01  1.98E+00  5.50E-01  1.51E+00  1.00E-02  1.28E+01  8.29E+00
 


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
+        3.80E+06
 
 TH 2
+       -6.98E+04  6.70E+07
 
 TH 3
+       -3.85E+03  3.63E+07  1.10E+08
 
 TH 4
+        6.32E+03 -3.77E+04 -7.78E+05  9.14E+06
 
 TH 5
+       -7.92E+00  1.26E+02  1.42E+02 -9.33E+01  1.08E+01
 
 TH 6
+       -3.90E+03  2.69E+06  1.54E+03 -1.10E+03  6.81E+00  4.01E+05
 
 TH 7
+       -4.86E+01  1.50E+01  6.59E+00 -1.19E+00  1.10E-02  2.22E+00 -2.35E-01
 
 TH 8
+       -2.64E+03  3.82E+04  1.43E+04 -6.92E+03  3.26E+01 -5.93E+05  3.24E-02  1.65E+06
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        2.15E+01 -4.08E+02 -2.26E+03  1.41E+03 -9.93E-01 -4.91E+01 -2.75E-02 -8.95E+01  0.00E+00  2.43E+02
 
 TH11
+       -6.34E+01  9.20E+02  2.51E+03 -1.27E+03  2.65E+00 -2.48E+04 -3.93E-03  2.24E+02  0.00E+00 -1.55E+01  8.13E+02
 
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
 #CPUT: Total CPU Time in Seconds,       45.806
Stop Time:
Sat Sep 25 14:29:22 CDT 2021

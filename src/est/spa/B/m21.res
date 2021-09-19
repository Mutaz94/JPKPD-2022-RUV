Sat Sep 18 08:21:35 CDT 2021
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
$DATA ../../../../data/spa/B/dat21.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1743.31078401763        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -6.2999E+01 -4.0020E+01 -9.9068E+00 -7.1393E+01 -4.9560E+01  1.8547E+00  3.5680E-01  1.6811E+01  2.8733E+00  2.0341E+01
             3.6271E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1753.60882123661        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0531E+00  1.0536E+00  1.1932E+00  1.0036E+00  1.1544E+00  9.7569E-01  1.0082E+00  8.7709E-01  9.7521E-01  9.3861E-01
             8.9697E-01
 PARAMETER:  1.5170E-01  1.5225E-01  2.7666E-01  1.0357E-01  2.4361E-01  7.5392E-02  1.0817E-01 -3.1142E-02  7.4894E-02  3.6641E-02
            -8.7357E-03
 GRADIENT:   9.6239E+01 -9.3130E+00  9.6679E+00 -2.6727E+01  1.1675E+01 -2.4968E+00 -1.8776E+00 -5.1746E-01 -6.5116E+00 -1.7891E+01
            -1.8948E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1754.78584598666        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0451E+00  1.0026E+00  1.0848E+00  1.0351E+00  1.0733E+00  9.9282E-01  1.0048E+00  6.0307E-01  9.6495E-01  9.4626E-01
             8.9826E-01
 PARAMETER:  1.4414E-01  1.0264E-01  1.8137E-01  1.3450E-01  1.7070E-01  9.2790E-02  1.0479E-01 -4.0573E-01  6.4320E-02  4.4765E-02
            -7.2998E-03
 GRADIENT:   7.3218E+01 -7.2775E+00  1.1405E+01 -1.6066E+01 -2.7827E+00  4.3214E+00 -5.1959E+00 -1.5591E+00 -4.4618E+00 -9.6253E+00
            -1.6796E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1755.19809941496        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0335E+00  1.0130E+00  9.8306E-01  1.0288E+00  1.0349E+00  9.9168E-01  1.1144E+00  5.0275E-01  9.2891E-01  9.3850E-01
             9.0551E-01
 PARAMETER:  1.3295E-01  1.1292E-01  8.2918E-02  1.2837E-01  1.3432E-01  9.1646E-02  2.0827E-01 -5.8767E-01  2.6260E-02  3.6525E-02
             7.4056E-04
 GRADIENT:   3.7846E+01 -7.2311E+00 -3.9737E+00 -5.0998E+00  6.1991E+00  3.0680E+00 -1.7679E+00  6.9079E-01 -2.3795E+00 -2.1004E+00
            -8.9141E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1755.98062237439        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      408
 NPARAMETR:  1.0520E+00  8.7110E-01  1.0823E+00  1.1385E+00  1.0144E+00  9.9784E-01  1.2726E+00  5.2241E-01  8.8085E-01  9.6194E-01
             9.2519E-01
 PARAMETER:  1.5065E-01 -3.7999E-02  1.7908E-01  2.2968E-01  1.1425E-01  9.7841E-02  3.4108E-01 -5.4931E-01 -2.6873E-02  6.1199E-02
             2.2240E-02
 GRADIENT:   6.6790E+00  6.4272E+00  3.1421E+00  7.8368E+00 -2.9364E+00  6.5176E-01  2.5515E-01 -5.4048E-01 -6.6956E-01 -7.8324E-01
            -8.0196E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1756.17280571013        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.0455E+00  6.9949E-01  1.1985E+00  1.2436E+00  9.9977E-01  9.9197E-01  1.4233E+00  6.7545E-01  8.4249E-01  9.8202E-01
             9.2766E-01
 PARAMETER:  1.4446E-01 -2.5740E-01  2.8108E-01  3.1801E-01  9.9772E-02  9.1935E-02  4.5299E-01 -2.9238E-01 -7.1392E-02  8.1852E-02
             2.4909E-02
 GRADIENT:  -1.7056E+00  1.2531E+00  1.1286E+00 -8.5617E-01 -4.8106E+00 -5.8036E-01 -1.4683E-01  4.7817E-01 -3.1344E-01  1.1467E+00
             1.7206E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1756.26816213639        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      760
 NPARAMETR:  1.0438E+00  5.6714E-01  1.3134E+00  1.3290E+00  1.0069E+00  9.9011E-01  1.5460E+00  7.7057E-01  8.2088E-01  9.9982E-01
             9.2759E-01
 PARAMETER:  1.4286E-01 -4.6715E-01  3.7259E-01  3.8440E-01  1.0692E-01  9.0056E-02  5.3568E-01 -1.6062E-01 -9.7381E-02  9.9815E-02
             2.4833E-02
 GRADIENT:  -2.4754E-01 -1.1237E-01  4.4400E-01 -2.1686E+00  1.0033E+00 -4.5538E-01  7.4734E-02 -2.4799E-01 -5.8847E-01 -4.8164E-01
             1.0436E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1756.45092976495        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      935
 NPARAMETR:  1.0398E+00  3.5575E-01  1.4352E+00  1.4688E+00  9.8353E-01  9.8834E-01  1.5685E+00  8.9120E-01  7.9824E-01  1.0221E+00
             9.1929E-01
 PARAMETER:  1.3903E-01 -9.3354E-01  4.6128E-01  4.8443E-01  8.3392E-02  8.8274E-02  5.5009E-01 -1.5187E-02 -1.2534E-01  1.2184E-01
             1.5851E-02
 GRADIENT:  -8.2154E-01  2.7899E+00  1.8489E+00  1.2590E+01 -5.4969E+00 -9.5760E-03 -1.4455E-04 -2.2387E-01 -9.5793E-01  6.9170E-01
            -2.1591E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1756.67866789326        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1114
 NPARAMETR:  1.0373E+00  1.9524E-01  1.6088E+00  1.5721E+00  1.0026E+00  9.8570E-01  1.2815E+00  1.0617E+00  7.6958E-01  1.0369E+00
             9.2572E-01
 PARAMETER:  1.3661E-01 -1.5335E+00  5.7549E-01  5.5242E-01  1.0258E-01  8.5593E-02  3.4804E-01  1.5984E-01 -1.6191E-01  1.3621E-01
             2.2814E-02
 GRADIENT:  -5.4315E-02  1.1613E+00  1.2547E+00  1.1265E+01 -1.3308E+00 -2.1889E-02  1.5048E-01 -2.7685E-01 -8.8022E-01 -5.4195E-01
             3.8110E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1756.76349285580        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1289
 NPARAMETR:  1.0360E+00  1.2849E-01  1.5842E+00  1.6050E+00  9.7669E-01  9.8468E-01  9.6525E-01  1.0398E+00  7.5669E-01  1.0325E+00
             9.2362E-01
 PARAMETER:  1.3541E-01 -1.9519E+00  5.6008E-01  5.7314E-01  7.6419E-02  8.4566E-02  6.4633E-02  1.3907E-01 -1.7880E-01  1.3202E-01
             2.0541E-02
 GRADIENT:   1.6907E-01 -2.9005E-02 -4.3221E-01 -6.5485E-01  2.7050E-01  1.0272E-02  6.1399E-02  2.4286E-02 -1.1165E-01  1.2476E-01
             3.5988E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1756.76897928768        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1464
 NPARAMETR:  1.0354E+00  1.0037E-01  1.6095E+00  1.6236E+00  9.7714E-01  9.8419E-01  8.0248E-01  1.0635E+00  7.4905E-01  1.0325E+00
             9.2283E-01
 PARAMETER:  1.3481E-01 -2.1989E+00  5.7593E-01  5.8464E-01  7.6878E-02  8.4067E-02 -1.2005E-01  1.6159E-01 -1.8895E-01  1.3196E-01
             1.9688E-02
 GRADIENT:  -3.9294E-02 -9.2852E-03  9.7372E-02  8.0003E-02 -8.0588E-02 -1.2965E-03  3.0005E-02 -2.3871E-02 -3.5612E-02 -2.4384E-02
            -3.9334E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1756.76941945363        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1649
 NPARAMETR:  1.0354E+00  1.0117E-01  1.6092E+00  1.6231E+00  9.7738E-01  9.8422E-01  7.9549E-01  1.0626E+00  7.4950E-01  1.0322E+00
             9.2278E-01
 PARAMETER:  1.3482E-01 -2.1909E+00  5.7575E-01  5.8432E-01  7.7125E-02  8.4090E-02 -1.2879E-01  1.6074E-01 -1.8835E-01  1.3170E-01
             1.9630E-02
 GRADIENT:  -4.7801E-02 -1.9718E-02  8.1937E-02  2.8598E-02  1.1595E-01  8.4697E-04  3.0144E-02 -7.1013E-02  6.4538E-02 -1.0675E-01
            -8.6375E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1756.79540076250        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1828
 NPARAMETR:  1.0365E+00  1.5471E-01  1.5972E+00  1.5923E+00  9.8945E-01  9.8498E-01  4.3026E-01  1.0559E+00  7.6380E-01  1.0349E+00
             9.2309E-01
 PARAMETER:  1.3590E-01 -1.7662E+00  5.6823E-01  5.6518E-01  8.9395E-02  8.4870E-02 -7.4336E-01  1.5440E-01 -1.6945E-01  1.3428E-01
             1.9976E-02
 GRADIENT:   6.7642E-02  2.7033E-01 -1.0023E+00  5.0052E+00  1.2472E+00 -6.2098E-02  2.0149E-02  5.3649E-02 -8.8898E-01 -4.4670E-01
            -2.3459E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1756.82570267816        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2004
 NPARAMETR:  1.0375E+00  2.0939E-01  1.5869E+00  1.5564E+00  1.0013E+00  9.8604E-01  2.1173E-01  1.0441E+00  7.8269E-01  1.0426E+00
             9.2387E-01
 PARAMETER:  1.3677E-01 -1.4636E+00  5.6179E-01  5.4240E-01  1.0135E-01  8.5938E-02 -1.4524E+00  1.4315E-01 -1.4502E-01  1.4171E-01
             2.0816E-02
 GRADIENT:  -1.5831E-01 -1.2688E-02 -1.5978E-01  1.7353E-01  3.5694E-01  1.0459E-02  9.9242E-03 -2.5318E-02 -1.6025E-01 -7.6148E-02
            -9.4320E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1756.83691921028        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2180
 NPARAMETR:  1.0385E+00  2.6488E-01  1.5756E+00  1.5221E+00  1.0138E+00  9.8689E-01  1.0223E-01  1.0345E+00  8.0125E-01  1.0482E+00
             9.2470E-01
 PARAMETER:  1.3781E-01 -1.2285E+00  5.5467E-01  5.2008E-01  1.1368E-01  8.6801E-02 -2.1805E+00  1.3396E-01 -1.2159E-01  1.4708E-01
             2.1713E-02
 GRADIENT:  -2.1885E-02  6.7254E-02  4.5285E-02  3.0979E-01 -1.2936E-01 -3.2745E-04  3.7012E-03 -1.0293E-02 -1.3332E-01 -1.0127E-02
            -2.6594E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1756.83811139132        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2368             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0388E+00  2.8040E-01  1.5740E+00  1.5125E+00  1.0179E+00  9.8713E-01  7.6702E-02  1.0329E+00  8.0672E-01  1.0503E+00
             9.2483E-01
 PARAMETER:  1.3810E-01 -1.1715E+00  5.5361E-01  5.1373E-01  1.1778E-01  8.7051E-02 -2.4678E+00  1.3241E-01 -1.1477E-01  1.4906E-01
             2.1850E-02
 GRADIENT:   7.3655E+01  4.6098E+00  1.0383E+00  1.0166E+02  5.7895E-01  5.5584E+00  7.4699E-03 -4.3406E-03  1.2692E+00  1.4476E-01
             4.1769E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1756.83897256956        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2545
 NPARAMETR:  1.0386E+00  2.7890E-01  1.5742E+00  1.5132E+00  1.0176E+00  9.8701E-01  3.6571E-02  1.0337E+00  8.0613E-01  1.0498E+00
             9.2487E-01
 PARAMETER:  1.3786E-01 -1.1769E+00  5.5376E-01  5.1421E-01  1.1748E-01  8.6928E-02 -3.2085E+00  1.3312E-01 -1.1551E-01  1.4861E-01
             2.1900E-02
 GRADIENT:  -4.3936E-01  6.7171E-03  4.4062E-02 -2.7019E-01 -3.6020E-02 -3.5365E-02  5.3150E-04 -5.1123E-03 -1.1248E-01 -3.1616E-02
            -8.9552E-04

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1756.83926137656        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2720
 NPARAMETR:  1.0387E+00  2.7849E-01  1.5725E+00  1.5134E+00  1.0170E+00  9.8707E-01  1.0000E-02  1.0321E+00  8.0629E-01  1.0496E+00
             9.2485E-01
 PARAMETER:  1.3795E-01 -1.1784E+00  5.5269E-01  5.1433E-01  1.1682E-01  8.6982E-02 -4.7769E+00  1.3157E-01 -1.1532E-01  1.4845E-01
             2.1879E-02
 GRADIENT:  -2.0753E-01 -1.0793E-02 -2.7316E-02 -2.4347E-01  4.1654E-03 -1.2243E-02  0.0000E+00  2.0373E-03 -1.0574E-03 -1.0254E-02
             1.4070E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1756.83926592184        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2897
 NPARAMETR:  1.0387E+00  2.7840E-01  1.5728E+00  1.5134E+00  1.0171E+00  9.8708E-01  1.0000E-02  1.0323E+00  8.0626E-01  1.0498E+00
             9.2487E-01
 PARAMETER:  1.3795E-01 -1.1787E+00  5.5288E-01  5.1438E-01  1.1691E-01  8.6991E-02 -5.2987E+00  1.3177E-01 -1.1535E-01  1.4858E-01
             2.1898E-02
 GRADIENT:  -2.0976E-01 -9.2211E-03 -2.5129E-02 -2.1980E-01  9.2452E-03 -8.9782E-03  0.0000E+00 -1.5308E-04  1.9992E-03  9.4357E-04
             2.0724E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1756.83928135590        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3080
 NPARAMETR:  1.0388E+00  2.7819E-01  1.5737E+00  1.5136E+00  1.0172E+00  9.8710E-01  1.0000E-02  1.0330E+00  8.0616E-01  1.0499E+00
             9.2480E-01
 PARAMETER:  1.3809E-01 -1.1794E+00  5.5341E-01  5.1452E-01  1.1709E-01  8.7011E-02 -6.8043E+00  1.3243E-01 -1.1548E-01  1.4867E-01
             2.1827E-02
 GRADIENT:   1.1763E-01  6.7578E-03  3.1680E-02 -1.1738E-01 -4.9985E-02  2.1528E-04  0.0000E+00 -8.8489E-03 -1.9914E-03  8.5914E-04
            -1.9979E-02

0ITERATION NO.:   99    OBJECTIVE VALUE:  -1756.83928418800        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     3211
 NPARAMETR:  1.0388E+00  2.7818E-01  1.5736E+00  1.5137E+00  1.0173E+00  9.8710E-01  1.0000E-02  1.0331E+00  8.0616E-01  1.0499E+00
             9.2482E-01
 PARAMETER:  1.3808E-01 -1.1795E+00  5.5340E-01  5.1453E-01  1.1710E-01  8.7012E-02 -6.8043E+00  1.3252E-01 -1.1548E-01  1.4866E-01
             2.1843E-02
 GRADIENT:   7.9678E-02  5.4874E-03  7.2326E-03 -1.0762E-01 -2.3079E-02  6.9668E-04  0.0000E+00 -2.9595E-03 -3.0986E-04  2.7240E-06
            -1.0413E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     3211
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5648E-04 -1.0772E-04 -2.8310E-02 -5.4782E-03 -3.7387E-02
 SE:             2.9857E-02  4.9903E-05  1.5231E-02  2.9374E-02  2.2037E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9582E-01  3.0890E-02  6.3075E-02  8.5205E-01  8.9780E-02

 ETASHRINKSD(%)  1.0000E-10  9.9833E+01  4.8973E+01  1.5932E+00  2.6173E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  7.3963E+01  3.1610E+00  4.5495E+01
 EBVSHRINKSD(%)  3.7236E-01  9.9843E+01  5.2328E+01  1.9712E+00  2.3159E+01
 EBVSHRINKVR(%)  7.4333E-01  1.0000E+02  7.7273E+01  3.9036E+00  4.0954E+01
 RELATIVEINF(%)  9.7526E+01  1.4024E-05  5.7171E+00  6.6557E+00  9.2607E+00
 EPSSHRINKSD(%)  4.4263E+01
 EPSSHRINKVR(%)  6.8934E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1756.8392841879963     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1021.6884576242581     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    41.67
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1756.839       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.78E-01  1.57E+00  1.51E+00  1.02E+00  9.87E-01  1.00E-02  1.03E+00  8.06E-01  1.05E+00  9.25E-01
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.10E+03
 
 TH 2
+       -2.07E+01  3.98E+02
 
 TH 3
+        1.18E+01  3.85E+01  8.01E+01
 
 TH 4
+        4.43E+01  4.76E+02 -3.11E+01  7.44E+02
 
 TH 5
+       -3.88E+01 -1.98E+02 -2.08E+02 -4.03E+01  6.00E+02
 
 TH 6
+        3.88E+01 -2.77E+01  1.32E+01 -4.13E+01 -1.71E+01  1.75E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -8.89E-01  6.28E-01 -5.69E+00  3.81E+00  5.07E+00 -4.25E+00  0.00E+00  2.69E+00
 
 TH 9
+        1.20E+02 -1.41E+02  1.82E+01 -3.53E+01 -1.33E+01  2.12E+01  0.00E+00 -3.56E+00  3.57E+02
 
 TH10
+        1.01E+01  1.72E+01  1.81E+01  3.51E+00 -6.35E+01 -6.38E+00  0.00E+00  2.41E+00  6.02E+00  1.04E+01
 
 TH11
+        3.56E+01  2.05E+00 -5.60E+01  5.33E+01  3.00E+01  1.65E+01  0.00E+00  3.10E+01  1.07E+01  3.07E+01  3.92E+02
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.05E+03
 
 TH 2
+       -2.46E+01  3.70E+02
 
 TH 3
+       -7.77E-02  4.64E+01  1.02E+02
 
 TH 4
+       -8.53E+00  4.51E+02 -2.11E+01  7.17E+02
 
 TH 5
+        3.42E+00 -1.97E+02 -1.91E+02 -3.85E+01  5.95E+02
 
 TH 6
+        1.29E+00 -2.15E+00 -3.61E-01 -2.30E+00 -3.01E+00  2.02E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -1.94E+00 -3.44E+00 -2.59E+01 -2.69E+00 -5.75E+00  4.52E-01  0.00E+00  2.67E+01
 
 TH 9
+        7.36E+00 -9.63E+01  7.52E+00 -8.32E-01 -1.43E+00  2.09E+00  0.00E+00 -5.49E-01  2.85E+02
 
 TH10
+        1.69E+00  8.71E+00 -2.08E+00 -2.68E-01 -6.80E+01  1.21E+00  0.00E+00  1.61E+01  4.26E+00  7.48E+01
 
 TH11
+       -7.84E+00 -1.42E+01 -1.77E+01 -1.03E+01 -8.48E+00  1.58E+00  0.00E+00  1.51E+01  2.29E+00  1.27E+01  2.56E+02
 
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
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.05E+03
 
 TH 2
+       -4.93E+01  3.48E+02
 
 TH 3
+       -1.59E+01  4.96E+01  9.99E+01
 
 TH 4
+       -6.74E+01  4.41E+02 -1.52E+01  7.17E+02
 
 TH 5
+        4.11E+01 -1.90E+02 -1.85E+02 -3.06E+01  5.94E+02
 
 TH 6
+       -3.27E+01  2.79E+01  6.33E+00  4.16E+01  1.61E+01  2.27E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        5.70E-01 -6.13E+00 -2.57E+01 -4.71E+00  2.20E+00 -9.77E+00  0.00E+00  1.98E+01
 
 TH 9
+       -8.24E+01 -5.89E+01  9.15E+00  2.77E+01  8.79E+00 -1.35E+01  0.00E+00 -5.62E+00  2.37E+02
 
 TH10
+        1.24E+01 -9.78E-01 -3.83E+00 -1.92E+01 -8.67E+01 -2.33E+01  0.00E+00  1.38E+01 -3.09E+00  8.26E+01
 
 TH11
+       -3.42E+01 -2.37E+01  6.25E+00 -5.05E+01 -3.25E+01 -5.54E+00  0.00E+00  4.80E+00  4.71E-01  1.57E+01  1.71E+02
 
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
 #CPUT: Total CPU Time in Seconds,       48.028
Stop Time:
Sat Sep 18 08:22:25 CDT 2021

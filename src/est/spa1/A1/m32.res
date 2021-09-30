Wed Sep 29 22:00:47 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat32.csv ignore=@
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
Current Date:       29 SEP 2021
Days until program expires : 200
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m32.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1687.02242309790        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4418E+02 -1.4527E+01  1.1629E+02 -5.4588E+01  6.8833E+01  7.6983E+01  5.3092E-01 -4.2054E+02 -9.3554E+01  3.5597E+00
            -3.8654E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1949.70623130729        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0043E+00  1.3661E+00  7.0863E-01  8.4405E-01  1.0466E+00  9.8627E-01  8.3592E-01  4.3652E+00  1.0039E+00  6.0271E-01
             1.1524E+00
 PARAMETER:  1.0427E-01  4.1199E-01 -2.4442E-01 -6.9540E-02  1.4553E-01  8.6180E-02 -7.9224E-02  1.5737E+00  1.0388E-01 -4.0633E-01
             2.4186E-01
 GRADIENT:   2.4829E+02  1.8676E+02 -2.8495E+01  7.9033E+01  2.1705E+01  6.1735E+01 -9.0462E+00  1.8726E+02 -1.3788E+01 -1.1838E+01
            -1.3343E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1964.68631916237        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      250
 NPARAMETR:  1.0108E+00  1.4669E+00  8.0016E-01  7.5032E-01  9.6543E-01  9.2915E-01  7.9257E-01  4.2933E+00  9.6581E-01  4.5634E-01
             1.2325E+00
 PARAMETER:  1.1079E-01  4.8313E-01 -1.2295E-01 -1.8726E-01  6.4815E-02  2.6516E-02 -1.3247E-01  1.5570E+00  6.5211E-02 -6.8451E-01
             3.0903E-01
 GRADIENT:  -8.0917E+01  8.4846E+00  1.3514E+01 -3.8233E+01 -1.2223E+02  1.3998E+01 -1.1124E+01  2.8698E+01 -3.3548E+01 -1.0813E+01
            -5.4178E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1969.21984911251        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      417
 NPARAMETR:  1.0109E+00  1.4656E+00  8.0002E-01  7.5032E-01  9.6560E-01  8.9205E-01  7.9257E-01  4.2692E+00  1.2405E+00  4.5685E-01
             1.2332E+00
 PARAMETER:  1.1080E-01  4.8228E-01 -1.2312E-01 -1.8725E-01  6.4999E-02 -1.4232E-02 -1.3247E-01  1.5514E+00  3.1553E-01 -6.8339E-01
             3.0965E-01
 GRADIENT:  -8.7298E+01  6.7341E+00  1.3480E+01 -1.9348E+01 -1.1178E+02 -1.5544E+00  1.0580E+01  2.8318E+01  2.4295E+00 -9.6979E+00
            -4.4869E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1971.63398301399        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  1.0107E+00  1.4648E+00  5.7723E-01  7.4979E-01  9.6630E-01  8.7110E-01  6.4996E-01  3.9793E+00  1.2586E+00  4.5700E-01
             1.2375E+00
 PARAMETER:  1.1060E-01  4.8169E-01 -4.4952E-01 -1.8796E-01  6.5721E-02 -3.7996E-02 -3.3084E-01  1.4811E+00  3.3000E-01 -6.8308E-01
             3.1313E-01
 GRADIENT:   2.0756E+02  2.2484E+02 -8.2962E+00  4.2329E+01 -2.0387E+01  1.2459E+01  5.0443E+00  1.7297E+02  2.9111E-01 -9.9311E+00
            -4.2823E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1973.20087935414        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  1.0106E+00  1.4647E+00  6.7297E-01  7.4977E-01  9.6630E-01  8.8975E-01  6.1435E-01  3.9731E+00  1.3444E+00  4.5699E-01
             1.2376E+00
 PARAMETER:  1.1059E-01  4.8163E-01 -2.9606E-01 -1.8799E-01  6.5714E-02 -1.6819E-02 -3.8719E-01  1.4795E+00  3.9595E-01 -6.8310E-01
             3.1316E-01
 GRADIENT:  -8.8625E+01 -8.7641E+00  1.1238E+00 -3.1914E+00 -6.0455E+01 -2.2796E+00  1.3797E-01  1.9872E+01 -4.4264E-01 -1.0879E+01
            -4.5070E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1974.08342509055        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      905             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0124E+00  1.4656E+00  6.2024E-01  7.4975E-01  9.6721E-01  8.7536E-01  5.7501E-01  3.6988E+00  1.3769E+00  4.6328E-01
             1.2451E+00
 PARAMETER:  1.1228E-01  4.8223E-01 -3.7764E-01 -1.8801E-01  6.6656E-02 -3.3116E-02 -4.5337E-01  1.4080E+00  4.1982E-01 -6.6943E-01
             3.1922E-01
 GRADIENT:   2.1453E+02  2.3395E+02  3.5886E-01  4.3877E+01 -2.2914E+01  1.4487E+01  5.3268E+00  1.5065E+02  1.4966E+01 -1.0888E+01
            -3.7673E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1976.24422244786        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1073
 NPARAMETR:  1.0403E+00  1.4676E+00  6.8579E-01  7.4892E-01  9.9447E-01  8.8477E-01  5.9549E-01  3.7122E+00  1.3784E+00  4.6248E-01
             1.2441E+00
 PARAMETER:  1.3950E-01  4.8365E-01 -2.7718E-01 -1.8913E-01  9.4450E-02 -2.2425E-02 -4.1837E-01  1.4116E+00  4.2095E-01 -6.7116E-01
             3.1842E-01
 GRADIENT:  -8.4839E+00 -3.0065E+01 -6.7258E-01 -1.0716E+00 -6.8650E+00 -1.9031E+00 -9.0075E-02 -1.5397E+00  2.5839E+00 -1.4523E+01
            -4.3148E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1976.70080230094        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1254
 NPARAMETR:  1.0431E+00  1.4686E+00  7.0243E-01  7.5055E-01  1.0019E+00  8.8620E-01  5.9839E-01  3.7126E+00  1.3778E+00  4.7229E-01
             1.2460E+00
 PARAMETER:  1.4220E-01  4.8432E-01 -2.5320E-01 -1.8696E-01  1.0187E-01 -2.0813E-02 -4.1351E-01  1.4117E+00  4.2048E-01 -6.5016E-01
             3.1995E-01
 GRADIENT:   1.6837E-01 -2.9825E+01  9.7775E-02  5.5633E-01 -3.3200E+00 -8.6170E-01  1.3967E-01 -3.4481E+00  3.0059E+00 -1.4875E+01
            -4.1233E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1976.82766181120        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1437
 NPARAMETR:  1.0390E+00  1.4770E+00  7.2281E-01  7.4720E-01  1.0123E+00  8.9774E-01  5.9254E-01  3.7513E+00  1.3807E+00  4.7137E-01
             1.2451E+00
 PARAMETER:  1.3824E-01  4.9002E-01 -2.2462E-01 -1.9143E-01  1.1222E-01 -7.8756E-03 -4.2334E-01  1.4221E+00  4.2256E-01 -6.5210E-01
             3.1924E-01
 GRADIENT:  -1.0642E+01 -2.6511E+01  1.3370E+00  1.6606E+00  7.7440E-01  4.0624E+00 -5.3201E-01 -4.2278E+00  2.1561E+00 -1.6237E+01
            -4.3273E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1980.44237932085        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1614
 NPARAMETR:  1.0292E+00  1.6872E+00  5.6989E-01  6.0950E-01  1.0722E+00  9.0950E-01  7.2156E-01  4.1120E+00  1.3664E+00  5.0761E-01
             1.3060E+00
 PARAMETER:  1.2873E-01  6.2308E-01 -4.6231E-01 -3.9511E-01  1.6976E-01  5.1425E-03 -2.2634E-01  1.5139E+00  4.1220E-01 -5.7803E-01
             3.6700E-01
 GRADIENT:  -3.9383E+01 -1.6533E+01 -4.5740E+00 -6.4993E+00 -1.5541E+01  7.9711E+00  5.1361E+00  5.0333E+00 -4.9918E+00 -1.9200E+01
            -9.5758E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1981.47658873156        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1797
 NPARAMETR:  1.0351E+00  1.7424E+00  5.5584E-01  5.7878E-01  1.0995E+00  8.8725E-01  7.1309E-01  4.1825E+00  1.3639E+00  5.1873E-01
             1.3243E+00
 PARAMETER:  1.3451E-01  6.5526E-01 -4.8728E-01 -4.4684E-01  1.9482E-01 -1.9632E-02 -2.3814E-01  1.5309E+00  4.1035E-01 -5.5638E-01
             3.8085E-01
 GRADIENT:  -2.6160E+01 -6.0240E+00 -4.2681E+00 -4.9194E+00 -9.1603E+00 -1.0457E+00  2.0169E+00  3.0987E+00 -8.7909E+00 -2.1702E+01
            -4.3793E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1981.85823030163        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1977
 NPARAMETR:  1.0365E+00  1.7522E+00  5.5790E-01  5.8074E-01  1.1003E+00  8.5982E-01  7.1331E-01  4.1247E+00  1.3704E+00  5.3261E-01
             1.3260E+00
 PARAMETER:  1.3589E-01  6.6086E-01 -4.8357E-01 -4.4345E-01  1.9561E-01 -5.1036E-02 -2.3783E-01  1.5170E+00  4.1511E-01 -5.2997E-01
             3.8219E-01
 GRADIENT:  -2.4164E+01  1.0757E+01 -3.1731E+00  2.2709E+00 -1.2931E+01 -1.3612E+01  2.2658E+00  2.9785E-01 -8.0044E+00 -2.1326E+01
            -1.8561E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1994.56933066914        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2139
 NPARAMETR:  1.0365E+00  1.7506E+00  5.5786E-01  5.8065E-01  1.1004E+00  8.8072E-01  6.9419E-01  4.1196E+00  1.4350E+00  8.0330E-01
             1.3261E+00
 PARAMETER:  1.3589E-01  6.5997E-01 -4.8365E-01 -4.4361E-01  1.9568E-01 -2.7017E-02 -2.6501E-01  1.5158E+00  4.6114E-01 -1.1902E-01
             3.8222E-01
 GRADIENT:  -2.2584E+01  4.0523E+00 -1.0047E+01  1.5626E+01 -6.6177E+01 -4.7225E+00 -9.9309E-02  3.5629E+01 -1.6486E+00 -9.8796E-01
             9.0207E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1994.59599461380        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:     2259
 NPARAMETR:  1.0365E+00  1.7503E+00  5.5787E-01  5.8064E-01  1.1004E+00  8.8604E-01  6.9299E-01  4.1185E+00  1.4359E+00  8.0471E-01
             1.3274E+00
 PARAMETER:  1.3589E-01  6.5977E-01 -4.8363E-01 -4.4363E-01  1.9569E-01 -2.0993E-02 -2.6674E-01  1.5155E+00  4.6183E-01 -1.1727E-01
             3.8325E-01
 GRADIENT:   3.1215E+02  4.1250E+02 -6.6776E+00  7.4518E+01 -5.8353E+01  2.0870E+01  7.1749E+00  7.8844E+01  9.9770E+00 -4.6876E-01
             1.3645E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1994.63699394387        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:     2379
 NPARAMETR:  1.0365E+00  1.7496E+00  5.5788E-01  5.8061E-01  1.1004E+00  8.8787E-01  6.9388E-01  4.1166E+00  1.4359E+00  8.0426E-01
             1.3249E+00
 PARAMETER:  1.3588E-01  6.5936E-01 -4.8362E-01 -4.4367E-01  1.9572E-01 -1.8929E-02 -2.6546E-01  1.5150E+00  4.6182E-01 -1.1784E-01
             3.8136E-01
 GRADIENT:   3.1353E+02  4.1263E+02 -6.6892E+00  7.4390E+01 -5.8008E+01  2.1754E+01  7.2821E+00  7.8974E+01  1.0100E+01 -5.8257E-01
             1.2361E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1994.67543367766        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:     2499
 NPARAMETR:  1.0365E+00  1.7492E+00  5.5788E-01  5.8059E-01  1.1005E+00  8.8998E-01  6.8892E-01  4.1157E+00  1.4359E+00  8.0326E-01
             1.3196E+00
 PARAMETER:  1.3588E-01  6.5918E-01 -4.8360E-01 -4.4371E-01  1.9575E-01 -1.6560E-02 -2.7263E-01  1.5148E+00  4.6183E-01 -1.1907E-01
             3.7734E-01
 GRADIENT:   3.1644E+02  4.1630E+02 -6.7091E+00  7.4579E+01 -5.7657E+01  2.2851E+01  6.7061E+00  7.9298E+01  9.4989E+00 -8.5059E-01
             9.2827E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1994.74646897487        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     2643
 NPARAMETR:  1.0365E+00  1.7488E+00  5.5791E-01  5.8054E-01  1.1006E+00  8.9052E-01  6.9428E-01  4.1079E+00  1.4360E+00  8.0348E-01
             1.3207E+00
 PARAMETER:  1.3588E-01  6.5894E-01 -4.8355E-01 -4.4380E-01  1.9581E-01 -1.5947E-02 -2.6488E-01  1.5129E+00  4.6183E-01 -1.1880E-01
             3.7814E-01
 GRADIENT:   3.1595E+02  4.1397E+02 -6.6516E+00  7.4349E+01 -5.7288E+01  2.3046E+01  7.3262E+00  7.8934E+01  1.0159E+01 -8.0840E-01
             1.0083E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1994.74797578477        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:     2745
 NPARAMETR:  1.0365E+00  1.7483E+00  5.5792E-01  5.8053E-01  1.1006E+00  8.8671E-01  6.8910E-01  4.1067E+00  1.4360E+00  8.0353E-01
             1.3209E+00
 PARAMETER:  1.3587E-01  6.5866E-01 -4.8355E-01 -4.4382E-01  1.9582E-01 -2.0235E-02 -2.7237E-01  1.5126E+00  4.6183E-01 -1.1874E-01
             3.7834E-01
 GRADIENT:  -2.2157E+01  1.5122E+00 -1.0142E+01  1.4321E+01 -6.4758E+01 -2.0596E+00 -7.2660E-01  3.5297E+01 -2.2441E+00 -1.1091E+00
             6.2763E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1994.83994940199        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     2884
 NPARAMETR:  1.0365E+00  1.7479E+00  5.5795E-01  5.8048E-01  1.1006E+00  8.9117E-01  6.9650E-01  4.0998E+00  1.4360E+00  8.0252E-01
             1.3156E+00
 PARAMETER:  1.3587E-01  6.5844E-01 -4.8349E-01 -4.4391E-01  1.9590E-01 -1.5222E-02 -2.6168E-01  1.5109E+00  4.6184E-01 -1.1999E-01
             3.7429E-01
 GRADIENT:  -2.1736E+01  3.8585E-01 -1.0073E+01  1.4036E+01 -6.4466E+01 -1.4106E-01  2.0303E-01  3.5056E+01 -1.4393E+00 -1.4079E+00
             3.6530E+00

0ITERATION NO.:   96    OBJECTIVE VALUE:  -1994.83994940199        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     2908
 NPARAMETR:  1.0365E+00  1.7479E+00  5.5795E-01  5.8048E-01  1.1006E+00  8.9117E-01  6.9650E-01  4.0998E+00  1.4360E+00  8.0252E-01
             1.3156E+00
 PARAMETER:  1.3587E-01  6.5844E-01 -4.8349E-01 -4.4391E-01  1.9590E-01 -1.5222E-02 -2.6168E-01  1.5109E+00  4.6184E-01 -1.1999E-01
             3.7429E-01
 GRADIENT:  -2.2631E+01 -5.4046E+00 -3.3675E+00  6.4353E+00 -5.7949E+01 -1.2650E-01  7.0050E+02 -2.2374E-01 -1.4779E+00 -1.5266E+03
            -9.6535E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2908
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0698E-02 -4.2476E-02 -2.1900E-02  2.3902E-02 -4.2620E-02
 SE:             2.9803E-02  2.0838E-02  1.7192E-02  2.3473E-02  1.7427E-02
 N:                     100         100         100         100         100

 P VAL.:         7.1963E-01  4.1512E-02  2.0272E-01  3.0854E-01  1.4458E-02

 ETASHRINKSD(%)  1.5568E-01  3.0190E+01  4.2404E+01  2.1364E+01  4.1618E+01
 ETASHRINKVR(%)  3.1112E-01  5.1265E+01  6.6827E+01  3.8163E+01  6.5916E+01
 EBVSHRINKSD(%)  6.7776E-01  2.8329E+01  2.7203E+01  2.4634E+01  4.2914E+01
 EBVSHRINKVR(%)  1.3509E+00  4.8632E+01  4.7006E+01  4.3199E+01  6.7412E+01
 RELATIVEINF(%)  9.8601E+01  9.0813E+00  2.7286E+01  1.1075E+01  1.3109E+01
 EPSSHRINKSD(%)  3.5815E+01
 EPSSHRINKVR(%)  5.8803E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1994.8399494019925     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1075.9014161973198     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    51.73
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.15
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1994.840       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.75E+00  5.58E-01  5.80E-01  1.10E+00  8.91E-01  6.97E-01  4.10E+00  1.44E+00  8.03E-01  1.32E+00
 


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
+        9.24E+05
 
 TH 2
+       -1.21E+01  4.02E+02
 
 TH 3
+        2.40E+05 -2.96E+04  2.52E+05
 
 TH 4
+       -1.03E+01  3.14E+04 -9.14E+01  2.77E+05
 
 TH 5
+       -3.02E+05  3.71E+04 -1.57E+05  1.10E+02  3.94E+05
 
 TH 6
+        2.49E+01 -2.36E+00  1.36E-01 -2.33E+00  6.48E-02  2.45E+02
 
 TH 7
+        1.65E+00 -2.70E+01  5.78E+00  6.79E+00 -1.35E+01  1.87E+02  2.76E+05
 
 TH 8
+        2.79E-01 -3.20E+00 -8.87E+00  5.38E+03 -1.78E+01  1.34E-01 -5.35E-01  4.83E+02
 
 TH 9
+        1.02E+00 -7.47E+00 -5.18E-01  5.33E+04  1.18E+00 -3.03E-01  4.11E+01  2.08E+03  4.17E+04
 
 TH10
+       -6.75E+05 -1.54E+01 -2.02E+00  1.27E+01 -5.41E+01  1.07E+06 -1.13E+03  1.31E+00  3.40E+00  9.87E+05
 
 TH11
+       -1.31E+05 -1.69E+01  9.09E-01  1.71E+00  1.16E+01 -6.70E+01 -2.07E+02  4.18E-01  6.86E+00 -3.62E+03  3.72E+04
 
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
 #CPUT: Total CPU Time in Seconds,       60.974
Stop Time:
Wed Sep 29 22:01:50 CDT 2021

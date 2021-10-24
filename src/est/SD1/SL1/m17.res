Sat Oct 23 14:45:09 CDT 2021
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
$DATA ../../../../data/SD1/SL1/dat17.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
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
 (2E4.0,E20.0,E4.0,2E2.0)

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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2811.76508210669        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6490E+02 -6.1781E+00  2.0458E+02  5.8320E+01  5.0948E+01  4.2808E+01 -7.1954E+01 -1.5190E+02 -4.1775E+01 -5.2887E+01
            -1.8777E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3262.93146789662        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0107E+00  1.0897E+00  7.8884E-01  9.5610E-01  9.8240E-01  9.9240E-01  1.1867E+00  1.1256E+00  1.0128E+00  1.2198E+00
             1.8094E+00
 PARAMETER:  1.1065E-01  1.8591E-01 -1.3720E-01  5.5110E-02  8.2238E-02  9.2370E-02  2.7115E-01  2.1835E-01  1.1269E-01  2.9865E-01
             6.9299E-01
 GRADIENT:   8.0444E+01 -5.3441E+00 -2.1718E+01 -3.4786E+01 -1.7164E+01  6.2256E+00  1.0654E+01  6.5023E-01  5.7584E+00  9.6789E+00
             1.3827E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3269.40669683067        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      262
 NPARAMETR:  1.0631E+00  1.5806E+00  1.3389E+00  7.0193E-01  1.7510E+00  9.5039E-01  9.1403E-01  2.0575E+00  1.1431E+00  1.9098E+00
             1.8195E+00
 PARAMETER:  1.6116E-01  5.5781E-01  3.9187E-01 -2.5393E-01  6.6016E-01  4.9112E-02  1.0104E-02  8.2149E-01  2.3377E-01  7.4699E-01
             6.9855E-01
 GRADIENT:   4.5437E+01 -8.8832E+01  2.4471E+01 -1.9108E+00  5.2105E+01 -2.6367E+01  1.3174E+00 -2.2924E+01  3.0887E+00  4.6186E+01
             1.4840E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3276.23025379627        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      451             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0486E+00  1.5928E+00  1.3615E+00  6.9165E-01  1.7853E+00  9.9753E-01  8.9188E-01  2.2351E+00  1.1134E+00  1.6119E+00
             1.8101E+00
 PARAMETER:  1.4745E-01  5.6552E-01  4.0862E-01 -2.6868E-01  6.7959E-01  9.7530E-02 -1.4418E-02  9.0430E-01  2.0738E-01  5.7742E-01
             6.9339E-01
 GRADIENT:   2.1942E+02  1.5644E+02  1.6151E+01  3.5268E+01  1.5227E+02  9.4623E+00 -4.5070E+00 -1.7051E+01  2.4785E+00  1.9766E+01
             1.5269E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3276.62223710839        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      587
 NPARAMETR:  1.0430E+00  1.5929E+00  1.3513E+00  6.9163E-01  1.7853E+00  1.0212E+00  9.4776E-01  2.2350E+00  1.0789E+00  1.6043E+00
             1.8100E+00
 PARAMETER:  1.4211E-01  5.6553E-01  4.0108E-01 -2.6870E-01  6.7958E-01  1.2101E-01  4.6345E-02  9.0424E-01  1.7595E-01  5.7272E-01
             6.9331E-01
 GRADIENT:  -1.0534E-01 -9.1708E+01  1.3483E+01  3.7853E-01  8.0538E+01  3.2016E+00  7.3475E-01 -2.0221E+01  2.9744E-01  1.0777E+00
             1.4079E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3279.51357929353        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      765
 NPARAMETR:  1.0453E+00  1.5949E+00  1.0042E+00  6.9119E-01  1.7825E+00  1.0148E+00  9.5300E-01  2.2332E+00  1.0099E+00  1.6981E+00
             1.8026E+00
 PARAMETER:  1.4432E-01  5.6680E-01  1.0423E-01 -2.6935E-01  6.7803E-01  1.1472E-01  5.1855E-02  9.0344E-01  1.0990E-01  6.2952E-01
             6.8926E-01
 GRADIENT:   3.8869E+00 -1.0028E+02 -2.5044E+00  2.4681E+01  1.0197E+02  6.1617E-01 -5.2349E-01 -2.6088E+00  8.6966E-02  6.6098E-01
             1.4776E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3285.87842663968        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      947
 NPARAMETR:  1.0447E+00  1.5974E+00  1.0024E+00  6.9162E-01  1.7790E+00  1.0138E+00  9.5588E-01  2.4943E+00  1.0081E+00  1.7228E+00
             1.6656E+00
 PARAMETER:  1.4377E-01  5.6838E-01  1.0236E-01 -2.6871E-01  6.7605E-01  1.1369E-01  5.4873E-02  1.0140E+00  1.0809E-01  6.4396E-01
             6.1021E-01
 GRADIENT:   5.9194E+00 -9.2831E+01 -3.4925E+00  2.6944E+01  9.6078E+01 -5.0202E-01  1.2146E-01  2.4024E-01  2.9095E-01 -5.3468E-01
            -1.8861E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3294.56152324166        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1127
 NPARAMETR:  1.0492E+00  1.6809E+00  8.7246E-01  6.9212E-01  1.5330E+00  1.0192E+00  9.3448E-01  2.5012E+00  1.0253E+00  1.6712E+00
             1.6629E+00
 PARAMETER:  1.4799E-01  6.1934E-01 -3.6442E-02 -2.6799E-01  5.2721E-01  1.1905E-01  3.2239E-02  1.0168E+00  1.2498E-01  6.1353E-01
             6.0856E-01
 GRADIENT:   1.4261E+01  1.9840E+01  5.9315E-01  5.0799E+01  7.9269E+00  1.2861E+00 -2.6121E-01  1.0744E+01  1.0729E+00  1.1133E+00
             1.6727E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3294.96410265138        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1309
 NPARAMETR:  1.0321E+00  1.6782E+00  8.6781E-01  6.8824E-01  1.5269E+00  1.0108E+00  9.4192E-01  2.4548E+00  9.8085E-01  1.6582E+00
             1.6709E+00
 PARAMETER:  1.3160E-01  6.1773E-01 -4.1786E-02 -2.7362E-01  5.2324E-01  1.1072E-01  4.0168E-02  9.9803E-01  8.0664E-02  6.0575E-01
             6.1338E-01
 GRADIENT:  -2.1146E+01  1.5818E+01  1.1714E+00  4.4950E+01  6.0239E+00 -2.0908E+00 -1.0802E+00  9.8015E+00 -2.0858E+00 -5.3715E-01
             9.4822E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3296.70210822262        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1470
 NPARAMETR:  1.0439E+00  1.6738E+00  8.4956E-01  6.5321E-01  1.5141E+00  1.0122E+00  9.3610E-01  2.3351E+00  1.0340E+00  1.6490E+00
             1.6658E+00
 PARAMETER:  1.4294E-01  6.1508E-01 -6.3039E-02 -3.2586E-01  5.1481E-01  1.1212E-01  3.3971E-02  9.4806E-01  1.3348E-01  6.0018E-01
             6.1033E-01
 GRADIENT:   2.3981E+02  3.1773E+02  4.1394E+00  4.7396E+01  6.0470E+01  1.9112E+01  5.7812E+00  1.3287E+01  2.6185E+00  2.1314E+01
             4.5344E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3298.31297765755        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1628
 NPARAMETR:  1.0227E+00  1.7282E+00  8.2807E-01  6.3318E-01  1.5612E+00  1.0110E+00  8.9574E-01  2.1520E+00  1.1019E+00  1.6290E+00
             1.6757E+00
 PARAMETER:  1.2243E-01  6.4705E-01 -8.8653E-02 -3.5700E-01  5.4546E-01  1.1097E-01 -1.0107E-02  8.6642E-01  1.9707E-01  5.8799E-01
             6.1621E-01
 GRADIENT:  -4.0243E+01 -1.3501E+01  2.1452E+00  1.2264E+01  4.2410E-01 -2.4808E+00 -1.7765E+00  8.0672E-01  4.5476E-01 -2.1930E+00
             1.5793E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3298.95446169229        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1810
 NPARAMETR:  1.0275E+00  1.7433E+00  7.8197E-01  6.1042E-01  1.5737E+00  1.0169E+00  8.9766E-01  2.0855E+00  1.1246E+00  1.6473E+00
             1.6774E+00
 PARAMETER:  1.2710E-01  6.5580E-01 -1.4594E-01 -3.9360E-01  5.5345E-01  1.1671E-01 -7.9638E-03  8.3500E-01  2.1739E-01  5.9913E-01
             6.1725E-01
 GRADIENT:  -2.9770E+01 -3.2831E+01 -1.0327E+00  1.3945E+00  3.1481E-01  9.0487E-02  1.7079E-01 -1.9963E-01  1.8199E+00 -3.1951E-01
             2.5459E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3299.51971619854        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1986
 NPARAMETR:  1.0359E+00  1.7690E+00  7.6346E-01  5.9050E-01  1.5825E+00  1.0073E+00  8.9709E-01  2.1438E+00  1.1067E+00  1.6471E+00
             1.6755E+00
 PARAMETER:  1.3530E-01  6.7044E-01 -1.6989E-01 -4.2678E-01  5.5900E-01  1.0725E-01 -8.6005E-03  8.6259E-01  2.0137E-01  5.9899E-01
             6.1611E-01
 GRADIENT:  -1.2394E+01 -3.4163E+01 -1.1158E+00 -5.1330E+00 -5.7753E+00 -3.1975E+00  2.9760E-01  1.1438E+00 -5.5405E-01 -2.4639E+00
            -1.0218E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3300.09527869701        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2164
 NPARAMETR:  1.0243E+00  1.8446E+00  6.9923E-01  5.3955E-01  1.6358E+00  1.0164E+00  8.6470E-01  2.0909E+00  1.1967E+00  1.6758E+00
             1.6787E+00
 PARAMETER:  1.2406E-01  7.1228E-01 -2.5778E-01 -5.1703E-01  5.9211E-01  1.1627E-01 -4.5376E-02  8.3761E-01  2.7954E-01  6.1630E-01
             6.1801E-01
 GRADIENT:  -3.6104E+01 -3.5850E+01 -2.1832E+00 -7.9554E+00 -5.5020E+00 -3.8105E-01 -1.0368E+00 -1.0718E+00  4.5715E-01 -1.0322E+00
             1.6492E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3300.38802473780        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2326
 NPARAMETR:  1.0240E+00  1.8857E+00  6.9880E-01  5.3059E-01  1.6381E+00  1.0144E+00  8.6220E-01  2.0951E+00  1.1959E+00  1.6814E+00
             1.6811E+00
 PARAMETER:  1.2376E-01  7.3430E-01 -2.5839E-01 -5.3377E-01  5.9353E-01  1.1427E-01 -4.8267E-02  8.3961E-01  2.7888E-01  6.1960E-01
             6.1947E-01
 GRADIENT:  -3.7137E+01 -1.5842E-01  1.3924E+00 -4.5919E-01 -1.2542E+01 -1.2650E+00 -3.1881E-01 -2.2347E+00 -1.3104E+00 -8.3390E-01
             1.0779E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -3300.50565668600        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     2486
 NPARAMETR:  1.0257E+00  1.8921E+00  6.9703E-01  5.2829E-01  1.6467E+00  1.0122E+00  8.6127E-01  2.1051E+00  1.1963E+00  1.6821E+00
             1.6797E+00
 PARAMETER:  1.2537E-01  7.3771E-01 -2.6093E-01 -5.3811E-01  5.9875E-01  1.1217E-01 -4.9351E-02  8.4435E-01  2.7923E-01  6.2004E-01
             6.1864E-01
 GRADIENT:   1.6936E+02  5.0326E+02  2.7103E+00  5.3157E+01  6.4408E+01  1.8975E+01  4.0441E+00 -6.8138E-01  1.9635E+00  1.9975E+01
             1.0164E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -3300.86719373771        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     2643
 NPARAMETR:  1.0368E+00  1.8848E+00  6.9081E-01  5.2653E-01  1.6636E+00  1.0168E+00  8.6181E-01  2.1729E+00  1.2137E+00  1.6880E+00
             1.6803E+00
 PARAMETER:  1.3617E-01  7.3384E-01 -2.6989E-01 -5.4144E-01  6.0899E-01  1.1671E-01 -4.8726E-02  8.7606E-01  2.9370E-01  6.2357E-01
             6.1895E-01
 GRADIENT:   2.1038E+02  4.8587E+02 -5.1895E-01  5.2895E+01  7.2098E+01  2.1495E+01  4.9077E+00  1.6653E+00  4.0856E+00  2.0741E+01
             1.2663E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -3301.01881097031        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2820
 NPARAMETR:  1.0409E+00  1.9117E+00  7.2120E-01  5.1394E-01  1.7020E+00  1.0137E+00  8.5083E-01  2.3102E+00  1.2293E+00  1.7022E+00
             1.6774E+00
 PARAMETER:  1.4010E-01  7.4801E-01 -2.2684E-01 -5.6565E-01  6.3179E-01  1.1364E-01 -6.1543E-02  9.3734E-01  3.0648E-01  6.3191E-01
             6.1724E-01
 GRADIENT:  -1.9876E+00 -3.6677E+00  3.2362E+00  4.5524E-01 -3.9325E-01 -7.7577E-01 -5.3791E-01 -2.2498E-02 -8.8254E-01  7.3798E-02
            -2.1523E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -3301.23669061713        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2995
 NPARAMETR:  1.0433E+00  1.9969E+00  6.3117E-01  4.6023E-01  1.7596E+00  1.0170E+00  8.3006E-01  2.2873E+00  1.3461E+00  1.7248E+00
             1.6825E+00
 PARAMETER:  1.4234E-01  7.9161E-01 -3.6018E-01 -6.7603E-01  6.6509E-01  1.1681E-01 -8.6253E-02  9.2737E-01  3.9724E-01  6.4513E-01
             6.2028E-01
 GRADIENT:   2.7992E+00  1.8636E+00 -5.9742E-01  1.6589E+00  4.3714E-01  3.4198E-01  3.5556E-01  2.0236E-01  9.0567E-02  8.2528E-02
             5.8406E-01

0ITERATION NO.:   91    OBJECTIVE VALUE:  -3301.23669061713        NO. OF FUNC. EVALS.:  27
 CUMULATIVE NO. OF FUNC. EVALS.:     3022
 NPARAMETR:  1.0429E+00  2.0001E+00  6.3345E-01  4.5946E-01  1.7585E+00  1.0166E+00  8.2924E-01  2.2916E+00  1.3457E+00  1.7252E+00
             1.6823E+00
 PARAMETER:  1.4234E-01  7.9161E-01 -3.6018E-01 -6.7603E-01  6.6509E-01  1.1681E-01 -8.6253E-02  9.2737E-01  3.9724E-01  6.4513E-01
             6.2028E-01
 GRADIENT:   7.9676E-01 -3.5377E+03 -4.9633E-01  9.0207E-01  1.0517E+00  1.5091E-01  2.9082E-01 -3.0679E+03  7.1385E-02 -1.5589E-01
             6.0353E-01
 NUMSIGDIG:         2.3         2.3         1.6         2.2         2.6         2.2         1.6         2.3         2.7         3.1
                    3.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3022
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1486E-04 -2.7228E-02 -2.4527E-02  3.4520E-02 -2.3573E-02
 SE:             2.9727E-02  2.5648E-02  1.4058E-02  1.9071E-02  2.6908E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9423E-01  2.8842E-01  8.1042E-02  7.0291E-02  3.8101E-01

 ETASHRINKSD(%)  4.1153E-01  1.4075E+01  5.2903E+01  3.6108E+01  9.8539E+00
 ETASHRINKVR(%)  8.2137E-01  2.6169E+01  7.7819E+01  5.9179E+01  1.8737E+01
 EBVSHRINKSD(%)  6.7745E-01  1.4112E+01  5.4312E+01  4.0690E+01  7.5036E+00
 EBVSHRINKVR(%)  1.3503E+00  2.6233E+01  7.9126E+01  6.4823E+01  1.4444E+01
 RELATIVEINF(%)  9.8641E+01  1.3605E+01  1.0829E+01  5.5418E+00  4.0914E+01
 EPSSHRINKSD(%)  1.9275E+01
 EPSSHRINKVR(%)  3.4835E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3301.2366906171260     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1647.1473308487152     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.49
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3301.237       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.00E+00  6.31E-01  4.60E-01  1.76E+00  1.02E+00  8.30E-01  2.29E+00  1.35E+00  1.72E+00  1.68E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,      233.233
Stop Time:
Sat Oct 23 14:45:43 CDT 2021

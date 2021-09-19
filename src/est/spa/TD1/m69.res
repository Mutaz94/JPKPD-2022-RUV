Sat Sep 18 14:13:09 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat69.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1679.68730742906        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.7511E+01  1.5175E+01 -5.8233E+01  9.4633E+01  9.2367E+01 -9.2907E+00 -6.7036E+00  8.7456E+00 -1.3429E+01 -6.5877E+00
             1.1898E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1687.29868098799        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      136
 NPARAMETR:  1.0284E+00  9.9711E-01  1.0782E+00  9.6157E-01  9.6832E-01  1.0297E+00  1.0234E+00  9.5323E-01  1.0629E+00  9.9171E-01
             9.5729E-01
 PARAMETER:  1.2801E-01  9.7110E-02  1.7533E-01  6.0812E-02  6.7805E-02  1.2931E-01  1.2317E-01  5.2098E-02  1.6103E-01  9.1678E-02
             5.6351E-02
 GRADIENT:  -2.1891E+00  4.7733E+00 -1.8620E+00  5.1530E+00 -1.2806E+00 -8.6097E-01 -3.7175E-01  4.2953E+00 -5.9611E-01 -4.0737E+00
            -5.8857E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1687.84400283574        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      314
 NPARAMETR:  1.0266E+00  9.3157E-01  1.1013E+00  1.0033E+00  9.5478E-01  1.0317E+00  1.0071E+00  7.9351E-01  1.0571E+00  1.0312E+00
             9.6472E-01
 PARAMETER:  1.2623E-01  2.9118E-02  1.9649E-01  1.0329E-01  5.3729E-02  1.3116E-01  1.0709E-01 -1.3129E-01  1.5551E-01  1.3070E-01
             6.4083E-02
 GRADIENT:  -5.0275E+00  4.2681E+00  1.5410E+00  6.4691E+00 -2.0802E+00 -3.7288E-02 -1.2273E+00  2.4491E-01  1.6882E+00 -2.1661E-01
            -2.0086E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1688.23949408518        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      490
 NPARAMETR:  1.0304E+00  9.1970E-01  9.3183E-01  9.9518E-01  8.8412E-01  1.0321E+00  1.1601E+00  5.4692E-01  1.0023E+00  9.5626E-01
             9.6872E-01
 PARAMETER:  1.2992E-01  1.6297E-02  2.9391E-02  9.5166E-02 -2.3167E-02  1.3161E-01  2.4854E-01 -5.0345E-01  1.0234E-01  5.5275E-02
             6.8223E-02
 GRADIENT:   3.2224E-01 -2.6816E+00 -2.4241E+00 -1.9556E+00  1.4292E+00 -2.4546E-01 -3.3782E-01  8.1769E-01  1.5838E-01  1.0115E+00
             5.0372E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1688.43407490401        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  1.0310E+00  9.9373E-01  8.2272E-01  9.4473E-01  8.6272E-01  1.0340E+00  1.1358E+00  2.6603E-01  1.0221E+00  9.2053E-01
             9.6782E-01
 PARAMETER:  1.3053E-01  9.3709E-02 -9.5142E-02  4.3140E-02 -4.7660E-02  1.3342E-01  2.2737E-01 -1.2242E+00  1.2185E-01  1.7198E-02
             6.7290E-02
 GRADIENT:  -8.9304E-01  2.1620E+00  1.5064E+00  1.1087E+00 -3.0228E+00  5.7173E-02  3.3593E-02  1.0521E-01 -4.7034E-01 -2.1173E-01
            -1.9078E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1688.50982846674        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      842
 NPARAMETR:  1.0320E+00  1.0976E+00  7.7890E-01  8.7673E-01  8.9095E-01  1.0341E+00  1.0529E+00  1.1226E-01  1.0818E+00  9.3093E-01
             9.6761E-01
 PARAMETER:  1.3148E-01  1.9312E-01 -1.4988E-01 -3.1560E-02 -1.5469E-02  1.3354E-01  1.5158E-01 -2.0869E+00  1.7860E-01  2.8432E-02
             6.7071E-02
 GRADIENT:  -2.0132E-02 -5.6999E-01 -1.6161E-01 -7.9172E-01  3.7984E-01 -1.4058E-02 -4.1853E-02  3.3903E-02 -5.5263E-02 -4.8222E-02
             3.2722E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1688.52680476797        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1018
 NPARAMETR:  1.0320E+00  1.0727E+00  7.7761E-01  8.9234E-01  8.7804E-01  1.0341E+00  1.0739E+00  3.0530E-02  1.0651E+00  9.2275E-01
             9.6769E-01
 PARAMETER:  1.3149E-01  1.7021E-01 -1.5152E-01 -1.3903E-02 -3.0064E-02  1.3351E-01  1.7129E-01 -3.3890E+00  1.6307E-01  1.9604E-02
             6.7161E-02
 GRADIENT:   1.9066E-02  2.9220E-01 -9.8182E-03  1.6316E-01 -4.3287E-01 -3.9783E-02 -5.9889E-02  2.4499E-03 -1.2371E-01 -4.3093E-03
             8.4179E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1688.52832935263        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1193
 NPARAMETR:  1.0320E+00  1.0786E+00  7.7775E-01  8.8863E-01  8.8121E-01  1.0342E+00  1.0694E+00  1.0000E-02  1.0695E+00  9.2515E-01
             9.6751E-01
 PARAMETER:  1.3148E-01  1.7567E-01 -1.5135E-01 -1.8078E-02 -2.6463E-02  1.3361E-01  1.6711E-01 -4.8196E+00  1.6721E-01  2.2202E-02
             6.6974E-02
 GRADIENT:  -2.0220E-03 -3.8344E-03  1.0602E-03 -4.7384E-03  1.4114E-03  9.3837E-05  4.2221E-04  0.0000E+00  2.4248E-04 -2.3200E-04
            -1.4457E-04

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1688.52832935263        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1215
 NPARAMETR:  1.0320E+00  1.0786E+00  7.7775E-01  8.8863E-01  8.8121E-01  1.0342E+00  1.0694E+00  1.0000E-02  1.0695E+00  9.2515E-01
             9.6751E-01
 PARAMETER:  1.3148E-01  1.7567E-01 -1.5135E-01 -1.8078E-02 -2.6463E-02  1.3361E-01  1.6711E-01 -4.8196E+00  1.6721E-01  2.2202E-02
             6.6974E-02
 GRADIENT:  -2.0220E-03 -3.8344E-03  1.0602E-03 -4.7384E-03  1.4114E-03  9.3837E-05  4.2221E-04  0.0000E+00  2.4248E-04 -2.3200E-04
            -1.4457E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1215
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.5916E-04 -1.3201E-02 -3.7870E-04  5.9196E-03 -1.9783E-02
 SE:             2.9856E-02  2.0518E-02  1.6354E-04  2.4817E-02  2.3846E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9575E-01  5.1996E-01  2.0581E-02  8.1147E-01  4.0674E-01

 ETASHRINKSD(%)  1.0000E-10  3.1262E+01  9.9452E+01  1.6860E+01  2.0114E+01
 ETASHRINKVR(%)  1.0000E-10  5.2751E+01  9.9997E+01  3.0878E+01  3.6183E+01
 EBVSHRINKSD(%)  3.8195E-01  3.0799E+01  9.9519E+01  1.7408E+01  1.8932E+01
 EBVSHRINKVR(%)  7.6244E-01  5.2113E+01  9.9998E+01  3.1785E+01  3.4280E+01
 RELATIVEINF(%)  9.9010E+01  2.5347E+00  3.3713E-04  4.4432E+00  7.6005E+00
 EPSSHRINKSD(%)  4.4453E+01
 EPSSHRINKVR(%)  6.9145E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1688.5283293526293     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -953.37750278889109     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.51
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1688.528       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.08E+00  7.78E-01  8.89E-01  8.81E-01  1.03E+00  1.07E+00  1.00E-02  1.07E+00  9.25E-01  9.68E-01
 


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
+        9.70E+02
 
 TH 2
+       -5.48E+00  4.20E+02
 
 TH 3
+        1.39E+01  1.93E+02  4.83E+02
 
 TH 4
+       -1.10E+01  3.53E+02 -2.20E+02  8.44E+02
 
 TH 5
+       -2.15E+00 -3.52E+02 -5.75E+02  2.67E+02  1.08E+03
 
 TH 6
+        5.50E-01 -1.12E+00  1.60E+00 -5.64E+00 -1.38E+00  1.84E+02
 
 TH 7
+        2.91E-01  2.06E+01 -9.19E-01 -8.01E+00 -1.06E+01 -2.79E-01  4.53E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.10E+00 -2.17E+01 -2.54E+01  3.54E+01  4.50E+00 -2.26E-01  2.14E+01  0.00E+00  8.80E+01
 
 TH10
+       -1.86E+00 -1.19E+01 -5.58E+01 -1.43E+01 -5.36E+01  1.31E+00  1.18E+01  0.00E+00  7.42E+00  1.03E+02
 
 TH11
+       -5.06E+00 -1.35E+01 -3.24E+01 -3.40E+00  7.69E+00  3.97E-01  5.62E+00  0.00E+00  6.35E+00  1.92E+01  2.28E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.669
Stop Time:
Sat Sep 18 14:13:32 CDT 2021

%Muldoon JJ, Chuang Y, Bagheri N, Leonard JN
%Macrophages employ quorum licensing to regulate collective activation

%This model represents a single cell, or a uniform population (homogeneous case)


function sim = Macrophage_model_homogeneous


%****Perturbations and simulation conditions****%


%plating cell density: -0.247 for high, 62.6 for low
taudensity = -0.247;

%IL-10 pre-treatment: 1 with, 0 without
dIL10 = 0;

%sTNFR pre-treatment: 1 with, 0 without
dsTNFR = 0;

%LPS treatment: 1 with, 0 without
dLPS = 1;

%Dose of LPS treatment: a value of 1 indicates 100ng/ml
%Other values scale proportionally with dose, e.g., 0.1 for 10ng/ml.
wLPSdose = 1;

%Time of BFA treatment (in hours); or Inf if no treatment
tauBFA = Inf;

%Time points (in hours)
T = 0:0.1:24;


%****Initial value for a homogeneous population****%


%For plating at high density
ivNtot = 0.1;

%For plating at low density
%ivNtot = 0.058;


%****Estimated parameters****%


%Values of fitted parameters
KP = [
    %Parameters from the cell-intrisic model (without intercellular communication)
    19.4          %P1
    0.219         %P2
    6.53          %P3
    7.27          %P4
    %Parameters from the cell-extrinsic model (with intercellular communication)
    1.99          %P5
    1.70          %P6
    0.747         %P7
    4.17          %P8
    0.614 / 1     %P9 this value is scaled inversely by the number of simulated cells (one)
    0.234         %P10
    3.96e5        %P11
    7.30];        %P12
    
%Parameters fitted using the dynamical model
parvals = num2cell(KP);
[   kinactivateIKK,  ... %P1
    wNFkBtxNFkB,     ... %P2
    wNFkBtxIkB,      ... %P3
    wIL10txFBD,      ... %P4
    wNFkBtxTNF,      ... %P5
    wstabilize,      ... %P6
    wdestabilize,    ... %P7
    ksecretion,      ... %P8
    kactivateTNFR,   ... %P9
    tau_s,           ... %P10
    tau_IKKK,        ... %P11
    tau_ds           ... %P12
    ] = parvals{:};


%****Other parameters****%


%Global translation
ktl = 15;

%Global transcription
ktx = 1;

%NF-kB basal transcription (0 -> NFkBm)
wbasaltxNFkB = 0.0001 * ivNtot / 0.1;

%Delayed effect of IL-10 pre-treatent on transcription
tauIL10_1 = 6;
tauIL10_2 = 18;

%Parameters for population growth (independent of the density at plating)
r1 = 1.99;
r2 = 0.0365;

%A parameter involved in post-transcriptional destabilizing regulation
dsm = 3;


%****Prior knowledge parameters (sources cited in SI)****%


%Parameters for TLR4
kTLR4synthesis = 0.0185; %synthesis:   0 -> TLR4
kdegTLR4       = 0.185;  %degradation: TLR4 -> 0, TLR4* -> 0
kdegLPS        = 0.058;  %degradation: LPS -> 0
kactivateTLR4  = 10;     %activation:  TLR4 -> TLR4*, via LPS

%Parameters for kinases
r28 = 5.88;    %activation:   IKKK -> IKKK*, via TRAF6* from TLR4*
r29 = 150;     %inactivation: IKKK* -> IKKK %r
r30 = 60000;   %activation:   IKK -> IKK*, via IKKK*

%Params for NFkB
k13 = 2.1;     %degradation:  IkBR -> 0
k16 = 15;      %translation:  IkBR -> IkBc
k19 = 1.08;    %import:       IkBc -> IkBn
k23 = 0.72;    %export:       IkBn -> IkBc
k22 = 324;     %import:       NFkBc -> NFkBn
k30 = 49.7;    %export:       NFkB-IkBn -> NFkB-IkBc
k33 = 4.2;     %degradation:  IkBc -> 0
k36 = 4.2;     %degradation:  IkBn -> 0
k45 = 1800;    %association:  NFkBc + IkBc -> NFkB-IkB
k48 = 1800;    %association:  NFkBn + IkBn -> NFkB-IkBn
k57 = 81;      %association:  IKK* + IkBc -> IKK-IkBc
k60 = 666;     %association:  NFkB-IkBc + IKK* -> NFkB-IKK-IkBc
k63 = 1800;    %association:  IKK-IkBc + NFkBc -> NFkB-IKK-IkBc
k78 = 432;     %degradation:  NFkB-IKK-IkBc -> IKK* + NFkBc
Ne  = 1.5;     %export:       NFkBn -> NFkBc

%TNFR
kTNFRsynthesis = 0.0102;  %synthesis:   0 -> TNFR
kdegTNFR       = 0.102;   %degradation: TNFR -> 0, TNFR* -> 0
kTNFR_IKKK     = 0.588;   %activation:  IKKK -> IKKK*, via TNFR*

%NFkB
kdegNFkBm      = 0.14;    %degradation: NFkBm -> 0
kdegNFkB       = 0.36;    %degradation: NFkBc -> 0, NFkBn -> 0

%mCherry
kdegmCherrym   = 2.1;     %degradation: mCherrym -> 0
kdegmCherry    = 0.5;     %degradation: mCherry -> 0
kmature        = 0.7;     %maturation:  mCherry -> mCherryf

%TNF
kdegTNFm       = 2.1;     %degradation: TNFm -> 0
kdegTNF        = 0.7;     %degradation: TNF -> 0

%Extracellular
kdegTNFpool    = 0.27;    %degradation: TNFpool -> 0


%****Initial values***%


%Number of intracellular variables
nvar = 25;

%Initialize array of initial values
IV = zeros(nvar, 1);

%Initial value of RNA for NF-kB
ivNFkBm = 0.0007 * ivNtot / 0.1;

%Populate the IV vector with the initial values of the state variables
IV(1:nvar, 1) = [ ...
    0.1,        ...  %1  TLR4
    0,          ...  %2  TLR4*
    0.1,        ...  %3  TNFR
    0,          ...  %4  TNFR*
    0.1,        ...  %5  IKKK
    0,          ...  %6  IKKK*
    0.1,        ...  %7  IKK
    0,          ...  %8  IKK*
    ivNFkBm,    ...  %9  NFkBm
    0,          ...  %10 NFkBc
    0,          ...  %11 NFkBn
    0,          ...  %12 IkBm
    0,          ...  %13 IkBc
    0,          ...  %14 IkBn
    ivNtot,     ...  %15 NFkB-IkBc
    0,          ...  %16 NFkB-IkBn
    0,          ...  %17 IKK-IkB
    0,          ...  %18 NFkB-IKK-IkB
    1,          ...  %19 Stabilizing
    0,          ...  %20 Destabilizing
    0,          ...  %21 mCherrym
    0,          ...  %22 mCherry
    0,          ...  %23 mCherryf
    0,          ...  %24 TNFm
    0];              %25 TNF

%Initial values for extracellular variables
ExC = [0;             %TNFpool
    wLPSdose * dLPS]; %LPS

%Stabilizing and destabilizing regulation are adjusted with IL-10 pre-treatment
if dIL10 == 1
    %Stabilizing regulation (set to zero)
    IV((c - 1) * nvar + 19, 1) = 0;
    %Destabilizing regulation (set to arbitrarily large value)
    IV((c - 1) * nvar + 20, 1) = 999;
end


%****ODEs****%


[~, sim] = ode23s(@(t, y, options)[ ...
      
    % Y1 TLR4
      kTLR4synthesis                   ...                                 %synthesis:    0 -> TLR4
    - kactivateTLR4  * y(1)  .* y(27)  ...                                 %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y(1)                                                %degradation:  TLR4 -> 0
    
    % Y2 TLR4*
      kactivateTLR4  * y(1)  .* y(27)  ...                                 %activation:   TLR4 -> TLR4*, via LPS
    - kdegTLR4       * y(2)                                                %degradation:  TLR4* -> 0
    
    % Y3 TNFR
      kTNFRsynthesis                   ...                                 %synthesis:    0 -> TNFR
    - kactivateTNFR  * y(3)  .* y(26) * (1 - dsTNFR) ...                   %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y(3)                                                %degradation:  TNFR -> 0
    
    % Y4 TNFR*
      kactivateTNFR  * y(3)  .* y(26) * (1 - dsTNFR) ...                   %activation:   TNFR -> TNFR*, via TNFpool
    - kdegTNFR       * y(4)                                                %degradation:  TNFR* -> 0
    
    % Y5 IKKK
    - r28/10         * y(2)  .* y(5)  ...                                  %activation:   IKKK -> IKKK*, via TLR4*
    - kTNFR_IKKK     * y(4)  .* y(5)  ...                                  %activation:   IKKK -> IKKK*, via TNFR*
    + r29            * y(6)                                                %inactivation: IKKK* -> IKKK
    
    % Y6 IKKK*
      r28/10         * y(2)  .* y(5)  ...                                  %activation:   IKKK -> IKKK*, via TLR4*
    + kTNFR_IKKK     * y(4)  .* y(5)  ...                                  %activation:   IKKK -> IKKK*, via TNFR*
    - r29            * y(6)                                                %inactivation: IKKK* -> IKKK
   
    % Y7 IKK
    - r30            * y(6)  .* y(7)  ...                                  %activation:   IKK -> IKK*, via IKKK*
    + kinactivateIKK * y(8)                                                %inactivation: IKK* -> IKK
    
    % Y8 IKK*
      r30            * y(6)  .* y(7)  ...                                  %activation:   IKK -> IKK*, via IKKK*
    - kinactivateIKK * y(8)           ...                                  %inactivation: IKK* -> IKK
    - k57            * y(8)  .* y(13) ...                                  %association:  IKK* + IkBc -> IKK-IkBc
    - k60            * y(8)  .* y(15) ...                                  %association:  IKK* + NFkB-IkBc -> NFkB-IkB-IKKc
    + k78            * y(18)                                               %degradation:  NFkB-IkB-IKKc -> IKK + NFkBc
    
    % Y9 NFkBm
      ktx * (wbasaltxNFkB + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * y(11))  ...
    ./                 (1 + dLPS * t.^3 .* exp(-t) / 6 .* wNFkBtxNFkB * (y(11) + wIL10txFBD * dIL10)) ... %transcription: 0 -> NFkBm
    - kdegNFkBm * y(9)                                                     %degradation:  NFkBm -> 0
    
    % Y10 NFkBc
      ktl      * y(9)           ...                                        %translation:  0 -> NFkBc
    - k22      * y(10)          ...                                        %import:       NFkBc -> NFkBn
    + Ne       * y(11)          ...                                        %export:       NFkBn -> NFkBc
    - k45      * y(10) .* y(13) ...                                        %association:  NFkBc + IkBc -> NFkB-IkBc
    - k63      * y(10) .* y(17) ...                                        %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    + k78      * y(18)          ...                                        %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    - kdegNFkB * y(10)                                                     %degradation:  NFkBc -> 0
    
    % Y11 NFkBn
      k22      * y(10)          ...                                        %import:       NFkBc -> NFkBn
    - Ne       * y(11)          ...                                        %export:       NFkBn -> NFkBc
    - k48      * y(11) .* y(14) ...                                        %association:  NFkBn + IkBn -> NFkB-IkBn
    - kdegNFkB * y(11)                                                     %degradation:  NFkBc -> 0
    
    % Y12 IkBm
      ktx * wNFkBtxIkB * y(11)  ...
    ./ (1 + wNFkBtxIkB * y(11)) ...                                        %transcription: 0 -> IkBm
    - k13 * y(12)                                                          %degradation:  IkBm -> 0
    
    % Y13 IkBc
      k16 * y(12)          ...                                             %translation:  IkBm -> IkBc
    - k19 * y(13)          ...                                             %import:       IkBc -> IkBn
    + k23 * y(14)          ...                                             %export:       IkBn -> IkBc
    - k33 * y(13)          ...                                             %degradation:  IkBc -> 0
    - k45 * y(10) .* y(13) ...                                             %association:  NFkBc + IkBc -> NFkB-IkBc
    - k57 * y(8)  .* y(13)                                                 %association:  IKK* + IkBc -> IKK-IkBc
    
    % Y14 IkBn
      k19 * y(13)          ...                                             %import:       IkBc -> IkBn
    - k23 * y(14)          ...                                             %export:       IkBn -> IkBc
    - k36 * y(14)          ...                                             %degradation:  IkBn -> 0
    - k48 * y(11) .* y(14)                                                 %association:  NFkBn + IkBn -> NFkB-IkBn
    
    % Y15 NFkB-IkBc
      k30 * y(16)          ...                                             %export:       NFkB-IkBn -> NFkB-IkBc
    + k45 * y(10) .* y(13) ...                                             %association:  NFkBc + IkBc -> NFkB-IkBc
    - k60 * y(8)  .* y(15)                                                 %association:  IKK* + NFkB-IkBc -> NFkB-IKK-IkBc
    
    % Y16 NFkB-IkBn
    - k30 * y(16)          ...                                             %export:       NFkB-IkBn -> NFkB-IkBc
    + k48 * y(11) .* y(14)                                                 %association:  NFkBn + IkBn -> NFkB-IkBn
    
    % Y17 IKK-IkB
      k57 * y(8)  .* y(13) ...                                             %association:  IKK* + IkBc -> IkB-IKKc
    - k63 * y(17) .* y(10)                                                 %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    
    % Y18 NFkB-IKK-IkB
      k60 * y(8)  .* y(15) ...                                             %association:  NFkB-IkBc + IKK* -> NFkB-IkB-IKKc
    + k63 * y(10) .* y(17) ...                                             %association:  IkB-IKKc + NFkBc -> NFkB-IkB-IKKc
    - k78 * y(18)                                                          %degradation:  NFkB-IkB-IKKc -> IKK* + NFkBc
    
    % Y19 Stabilizing
    - (1 - dIL10) * tau_s * y(19) .* 1 ./ (1 + tau_IKKK * max(0, y(6) - 0.0000392))
    
    % Y20 Destabilizing
      1 - (1 - dIL10) * (1 - 1 ./ (1 + exp(-(t - tau_ds))))
    
    % Y21 mCherry
      dLPS * (ktx * (wNFkBtxTNF * y(11)   ...
    ./          (1 + wNFkBtxTNF * y(11))) ...                              %transcription: 0 -> mCherrym, via MAPKs* and NFkB
    * (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) ...
    .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...                 %delayed effect of IL-10 pre-treatment, if applicable
    - kdegmCherrym * y(21)                                                 %degradation: mCherrym -> 0
    
    % Y22 mCherry
      ktl         * y(21) ...                                              %translation: 0 -> mCherry
    - kdegmCherry * y(22)                                                  %degradation: mCherry -> 0
    
    % Y23 mCherryf
      kmature     * y(22) ...                                              %maturation:  mCherry -> mCherryf
    - kdegmCherry * y(23)                                                  %degradation: mCherryf -> 0
    
    % Y24 Tnfm
      dLPS * (ktx * wNFkBtxTNF * (y(11))  ...
    ./         (1 + wNFkBtxTNF * (y(11))) ...                              %transcription: 0 -> TNFm
    .* (1 - dIL10 * (1 - (1 - max(t - tauIL10_1, 0) / (tauIL10_2 - tauIL10_1)) .* (1 - 1 ./ (1 + exp(-99 .* (t - tauIL10_2))))))) ...  %delayed effect of IL-10 pre-treatment, if applicable
    - kdegTNFm * y(24) .* (1 ./ (1 + wstabilize * y(19)) + (dsm - 1) * wdestabilize .* y(20) ./ (1 + wdestabilize .* y(20)))  %degradation: TNFm -> 0, via MAPKs* and TTP*
    
    % Y25 TNF
      ktl        * y(24) .* ((1 + wdestabilize * y(20)) ./ (1 + dsm * wdestabilize * y(20))) ... %translation: 0 -> TNF
    - ksecretion * y(25) .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) ...    %secretion: TNF -> TNFpool
    - kdegTNF    * y(25)                                                   %degradation: TNF -> 0
      
    %****Extracellular variables****%
        
    % Y26 TNFpool
      r1 ./ (1 + exp(-r2 * (t - taudensity))) ...                          %secretion: TNF -> TNFpool
    .* (1 - (1 + exp(-99 * (t - tauBFA))).^-1) * ksecretion .* y(25) ...
    - kdegTNFpool * y(26)                                                  %degradation: TNFpool -> 0
    
    % Y27 LPS
    - kdegLPS * y(27)                                                      %degradation: LPS -> 0
    
    ], T,...    %Simulate over the specified timecourse
    [IV; ExC]); %Initial values
    


end


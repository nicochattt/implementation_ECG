%-------------------------------------------------------------------------------
%  Gnx_Probleme2.m
%
% Problème 2
% G5C - BULTE, CEBE, DEVILLE, GUILLON, HOU, JALABERT
%
%  Vous ne devez pas modifier l'architecture' du programme
%  Vous ne devez pas modifier le code existant.
%  Vous devez ajouter vote code dans les deux functions analyse_temp(...) et
%  analyse_freq(...).
%  Vous pouvez ajouter d'autres' functions si nécessaire.
%-------------------------------------------------------------------------------

clc;
close all;
clear variables;


FreqECG_T = zeros(1,10);
FreqECG_F = zeros(1,10);

%-------------------------------------------------------------------
% Dans un premier temps, tester un seul fichier avec : for nf = 1:1
% Quand l'algorithme marche, mettre la boucle pour les 10 fichiers'
%-------------------------------------------------------------------
for nf = 1:10
	NameFile      = "./Gnx_Probleme2/Fichiers_ECG/10"+num2str(nf-1)+".wav";
	[sigSrc, Fe]  = audioread(NameFile);

	FreqMoy       = analyse_temp(nf, sigSrc, Fe);
	FreqECG_T(nf) = FreqMoy;
	FreqMoy       = analyse_freq(nf, sigSrc, Fe);
	FreqECG_F(nf) = FreqMoy;
end 


ValRef  = [74 67 73 70 74 83 59 70 62 92];

disp('Numero Fichier : 0   1   2   3   4   5   6   7   8   9 ');
disp(['Freq attendues : ' num2str(ValRef)]);
disp('--------------------------------------------------------');
disp(['Freq par Pics  : ' num2str(FreqECG_T)]);
disp(['ecart          : ' num2str(FreqECG_T-ValRef,'%+d  ')]);
disp(' ');
disp(['Freq par FFT   : ' num2str(FreqECG_F)]);
disp(['ecart          : ' num2str(FreqECG_F-ValRef,'%+d  ')]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function freqBpm = analyse_temp(nf, sigSrc, Fe)
%-----------------------------------------------------------------------
% Determiner la frequence ECG par l'analyse temporelle'
% Utilisation de la puissance et/ou de la correlation
% Toute valeur numerique doit etre definie par une variable avec un
% commentaire clair et explicite.
% Remarque : supprimer les plots ajoutés dans la version finale à rendre
%-----------------------------------------------------------------------
	SecPerMn= 60;		% secondes par minute
	Te      = 1.0/Fe;   % Fréquence d'échantillon
    
    seuil_nul = 0.5;                                % Pour déterminer le début du signal
    indiceEcg = find(sigSrc > seuil_nul);
    sigEcg = sigSrc(indiceEcg(1):indiceEcg(end));   % Là où il y a effectivement un signal Ecg
    sigEcg = sigEcg - mean(sigEcg);                 % On retire la composante continue
    N = length(sigEcg);                             % Nombre d'échantillon

  
    % Déterminer la durée de la fenetre en fonctions des valeurs
	% d'un signal ECG' d'un' être humain.

    % Nombre de fenêtres pour calculer la fréquence cardiaque
    NbrFen = 4;                                % Nombre de fenêtre
    Nb_Ech_Fen = floor(N/NbrFen);              % Nombre d'échantillon par fenetre

    TabFreq = zeros(1, NbrFen);  % Fréquence de chaque fenêtre
    VecFreq = zeros(1, N);   % Vecteur de la taille du signal


    % Détermination de la fréquence
    for k = 1:NbrFen
        % On détermine la fenêtre
        n1 = (k-1)*Nb_Ech_Fen + 1;   % Indice de début de fenêtre
        n2 = min(k*Nb_Ech_Fen, N);   % Indice de fin de fenêtre
        fenetre = sigEcg(n1:n2);   % Fenetre keme

        % auto correlation
        signal_autocorr = xcorr(fenetre);
        % Identification des pics
        [Hpics, localisation] = findpeaks(signal_autocorr,'MinPeakHeight', 0.1); % trouver les pics avec une hauteur minimale de la moitié de la hauteur maximale du spectre


        % Affichage des amplitudes correspondant aux pics
        [~, max_idx] = max(Hpics);
        Hpics_demi_signal = Hpics(max_idx + 1:end);
        localisation_demi_signal = localisation(max_idx + 1:end);
        [~, max_idx_demi_signal] = max(Hpics_demi_signal);

        freqFen = 60/(abs(localisation_demi_signal(max_idx_demi_signal) - localisation(max_idx)) * Te );

        freqFen = round(freqFen);
        disp("nous avons une freq cardiaque de = " + freqFen);

        TabFreq(k) = freqFen;  
        VecFreq(n1:n2) = freqFen;
   
    end

	moyBpm	= mean(TabFreq);	% Calculer la frequence moyenne
	freqBpm = round(moyBpm);
	Affiche_Resultat(nf+10, sigSrc, sigEcg, VecFreq, freqBpm, Te)

end

function freqBpm = analyse_freq(nf, sigSrc, Fe)
%-----------------------------------------------------------------------
% Determiner la frequence ECG par l'analyse frequentielle'
% Calcul de FFT par fenetrage (au minimum 4 fenetres)
% Toute valeur numerique doit etre definie par une variable avec un
% commentaire clair et explicite.
% Remarque : supprimer les plots ajoutés dans la version finale à rendre
%-----------------------------------------------------------------------
	SecPerMn= 60;		% secondes par minute
	Te      = 1.0/Fe;   %Periode d'échantillonnage
    deltaBpm=0.5; %On veut une précision de 0,5bpm
	deltaF  = deltaBpm/60;		% précision < 1 battement /mn (on prend 0,5 bpm)


	%sigSrc = signal d'origine'
	sigEcg	= sigSrc(sigSrc>0.1); %Signal à traiter
	NbrFen	= 10; % Calculer la frequence dans M (min 4) fenetres

	% Déterminer la durée de la fenetre en fonctions des valeurs
	% d'un signal ECG' d'un' être humain.
    %Une fréquence cardiaque varie en général entre 50 et 180 bpm, prenons
    %entre 40 et 200 pour avoir de la marge. On souhaite avoir au moins 3
    %périodes par fenêtre, on a donc:
    nbPeriodesFenetres = 3;
    bpmMin=40;
    dureeFenetre=(nbPeriodesFenetres*bpmMin)/SecPerMn;


    tailleFenetre = dureeFenetre*Fe;
    tailleDemiFenetre = round(tailleFenetre/2);

    

	TabFreq= zeros(1,NbrFen); % mémoriser les freq des fenetres
	VecFreq= zeros(1,length(sigEcg));%Vecteur de la taille du signal
    frequences = (0:deltaF:Fe-deltaF);
    pas = round((length(sigEcg)-tailleFenetre)/NbrFen);  %Le pas utilisé pour décaler le centre de la fenêtre dans le signal

	k = 1;
    for centreFen = tailleDemiFenetre:pas:length(sigEcg)-tailleDemiFenetre
        n1 = centreFen-tailleDemiFenetre+1; %n1 et n2 sont les indices de début et de fin de la fenêtre respectivement
        n2 = centreFen+tailleDemiFenetre;

        valMoy=mean(sigEcg(n1:n2)); %On recentre notre signal
        sigEcg2=sigEcg(n1:n2)-valMoy;


        fftSig = abs(fft(sigEcg2.^2, length(frequences))); % On calcule la Fft du signal, puis son autocorrélation
        [corrFft,lags]= xcorr(fftSig);

        lenCorr = length(corrFft); %taille du vecteur d'autocorrélation
        semiLenCorr=round(lenCorr/2); %taille de la demi-autocorrélation, que l'on va étudier
        semiCorr=corrFft(semiLenCorr:end); %On extrait la moitié de l'autocorrélation
        
        [pics,freqFen] = findpeaks(semiCorr); %On trouve les pics de l'autocorrélation avec findpeaks
        

        HzPerBpm=60;
		%freqFen = frequence en bpm de la k eme fenetre
		TabFreq(k)	= (freqFen(1)*deltaF)*HzPerBpm; %La fondamentale est le premier pic rencontré.
		VecFreq(n1:n2) = (freqFen(1)*deltaF)*HzPerBpm;%On convertit l'indice du pic en Hz puis les Hz en Bpm.
		k = k+1;
    end

	moyBpm	= median(TabFreq);	% Calculer la frequence moyenne
	freqBpm = round(moyBpm);
	Affiche_Resultat(nf+10, sigSrc, sigEcg, VecFreq, freqBpm, Te)
end

function Affiche_Resultat(nf, sigSrc, sigEcg, VecFreq, freqBpm, Te)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
	figure(nf);
	if (nf <= 10)
		Titre = "Fichier  10"+num2str(nf-1)+"  : Analyse Temporelle";
	else
		Titre = "Fichier  10"+num2str(nf-11)+"  : Analyse Fréquentielle";
	end
	ln1 = length(sigSrc);	ln2 = length(sigEcg);
	t1  = (0:ln1-1)*Te;		t2  = (0:ln2-1)*Te;
 	subplot(3,1,1);			plot(t1, sigSrc);
	title(Titre);
 	subplot(3,1,2);			plot(t2, sigEcg);
	subplot(3,1,3);			plot(t2, VecFreq);
	xlabel('temps (s)');
	ylabel('Rythme C. (bpm)');
	%ylim([freqBpm*0.8 freqBpm*1.2]);
	Label = "Rythme C. = "+num2str(freqBpm)+" bpm  ";
	yline(freqBpm, '-.r', Label);
end

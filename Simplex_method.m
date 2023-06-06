% Sunoptiki perigrafi twn metavlitwn

% vars = plithos metavlitwn tis antikeimenikis sunartisis pou thelw na megistopoihsw
% vars =
%  3

% maxz = oroi tou maxz
% maxz =
%  4   2   1

% per = sthles x_1, x_2, x_3, oroi twn periorismwn
% per =
%  1   2   3
%  2   4   1

% new_per = sthles x_4, x_5, kainouries metavlites pou tha prostethoun stous periorismous kai stin antikeimeniki sunartisi, tha einai monadiaios pinakas stin arxi
% new_per =
%  1   0
%  0   1

% b = sthlh b_0, times pou den prepei na ksepernane oi periorismoi
% b =
%  20
%  10

% S = pinakas pou exei ta x_1, x_2, x_3, x_4, x_5, b_0, me 2 grammes pou aforoun tin metavliti basis
% S =
%  1   2   3   1   0   20
%  2   4   1   0   1   10

% cj = grammi cj, kostos antikeimenikis sunartisis, diladi oi oroi tou maxz
% cj =
%  4   2   1   0   0   0

% S_CjZj = pinakas pou exei ta x_1, x_2, x_3, x_4, x_5, b_0, me 3 grammes, oi 2 prwtes aforoun tin metavliti basis kai h teleutaia grammi to cj, diladi o enwmenos pinakas S me cj
% S_CjZj =
%  1   2   3   1   0   20
%  2   4   1   0   1   10
%  4   2   1   0   0   0 

format shortg

vars = 3;				% 3 metavlites apofasewn

maxz = [4 2 1];			% mi midenikoi oroi tou maxz

per = [1 2 3; 2 4 1];	% x_1, x_2, x_3 (periorismoi xwris to b)
% Oi grammes einai to posoi periorismoi einai, edw exoume 2 periorismous ara 2 grammes

new_per = eye(size(per,1));		% Ftiakse monadiaio pinaka megethous 2x2

b = [20; 10];		% b_0

S = [per new_per b];	% Pinakas me olous tous periorismous

cj = zeros(1, size(S,2)); 		% Ftiakse grammi megethous 5 kai vale times 0
cj(1:vars) = maxz;				% Tha einai h grammi cj

base = vars+1:1:size(S,2)-1;	% Oi kainouries metavlites tha einai oi metavlites basis x_B

CjZj = cj; 						% Cj-Zj grammi, stin arxi einai idia me tin Cj afou h Zj exei timi 0 se kathe sthlh

S_CjZj = [S ; CjZj];

Zj = zeros(1, size(S,2));

% Kanw tis sthles tou pinaka S_CjZj na ginontai metavlites
Simplex = array2table(S_CjZj);
Simplex.Properties.VariableNames(1:size(S_CjZj,2)) = {'x_1','x_2','x_3','x_4','x_5','b_0'};

str = strings(1,2);

for j=1:2

	chr(j)=int2str(base(j));

	str(j)=convertCharsToStrings(chr(j));

	str(j)=strcat("x"+str(j));

end

X_base1 = convertStringsToChars(str(1));
X_base2 = convertStringsToChars(str(2));

Simplex.Properties.RowNames={X_base1, X_base2,'cj-zj'};
fprintf(" ----------------- Arxikos Simplex pinakas ----------------- \n\n");
disp(Simplex);

flag = true;

while flag
	
	if any(CjZj>0);			% An toulaxiston 1 stoixeio tis Cj-Zj grammhs einai thetiko, tote den exei vrethei h veltisti lusi
		
		fprintf('\nDen exei vrethei h veltisti lusi.\n')

		new_CjZj = CjZj(1:end-1);					% Oli h Zj-Cj grammi xwris thn b_0 sthlh
		
		[max_CjZj, pivot_col_index] = max(new_CjZj);			% Pairnw to max stoixeio kai sthlh tou max stoixeiou tis Zj-Cj grammhs

		fprintf('Max stoixeio stin cj-Zj grammi: %d. Epomenws, sthlh-odhgos: %d.\n', max_CjZj, pivot_col_index);
		fprintf('Metabliti basis pou tha ginei h antikathistwsa: x_%d\n', pivot_col_index);

		if all(S(:,pivot_col_index)<=0)
			
			disp('To b_0 den mporouse na upologistei se kamia grammi. H methodos Simplex tha stamatisei.');
			
			flag = false;
		
		else % an ola ta stoixeia ths odhgou sthlhs einai arnitika h 0, den mporei na upologistei kamia timi tou b_0
			
			b0 = S(:,end);						% Ta stoixeia ths b_0 gia tis grammes twn periorismwn kai tis grammhs Cj
			
			pivot_col = S(:,pivot_col_index);	% Bazw ta stoixeia tis sthlh-odhgou se pinaka gia na ta xrisimopoihsw ston upologismo tou theta

			for i=1:size(S,1)					% Ypologise tin sthlh theta
				
				if pivot_col(i)>0

					th(i) = b0(i)./pivot_col(i);

				else

					th(i)=inf;					% An to stoixeio ths sthlhs-odhgou einai arnitiko, tote aprosdioristo to apotelesma kai den to lambanoume up' opsi

				end

			end

			[MinTh, pivot_row_index] = min(th);			% Pairnw to min stoixeio kai grammh tou min stoixeiou tis theta stilis
			
			fprintf('Min stoixeio stin τheta sthlh: %d. Ara, grammh-odhgos: %d\n', MinTh, pivot_row_index);

			fprintf('H metavliti vasis pou tha antikatastathei: x_%d\n', base(pivot_row_index));

		end

		base(pivot_row_index) = pivot_col_index;

		pivot_element = S(pivot_row_index, pivot_col_index);	% *** Stoixeio - odhgos, c)ii) erwtima ***

		S(pivot_row_index,:) = S(pivot_row_index,:)./pivot_element;			% Grammoprakseis sthn odhgo-grammh
		
		for i=1:size(S,1)

			if i~=pivot_row_index

				S(i,:)=S(i,:)-S(i,pivot_col_index).*S(pivot_row_index,:);		% Grammoprakseis stis grammes ektos tis grammis-odhgou

			end

		end

		CjZj = CjZj - CjZj(pivot_col_index).*S(pivot_row_index,:);				% Grammoprakseis stin grammi tou Cj-Zj

		S_CjZj = [S ; CjZj];

		% Kanw tis sthles tou pinaka S_CjZj na ginontai metavlites
		TABLE=array2table(S_CjZj);
		TABLE.Properties.VariableNames(1:size(S_CjZj,2)) = {'x_1','x_2','x_3','x_4','x_5','b_0'};

		Solution = zeros(1,size(S,2));			% Kanw pinaka megethous 1x6 kai arxikopoiw tis times tou me 0 
		Solution(base) = S(:,end);				% Ston deikti tis base kai meta, vazw ta stoixeia ths b_0 pou antistoixizontai me ta x_1, x_2, x_3, x_4, x_5
		Solution(end) = sum(Solution.*cj);		% Sto teleutaio stoixeio tou pinaka Solution vazw to athroisma twn b_0
		
		% Kanw tis sthles tou pinaka Solution na ginontai metavlites
		Current_Solution = array2table(Solution);
		Current_Solution.Properties.VariableNames(1:size(Current_Solution,2)) = {'x_1','x_2','x_3','x_4','x_5','b_0'};
		
		Final_simplex=[S ; CjZj];

		TABLE2=array2table(Final_simplex);
		TABLE2.Properties.VariableNames(1:size(S_CjZj,2)) = {'x_1','x_2','x_3','x_4','x_5','b_0'};
		
		
		for j=1:2

				chr(j)=int2str(base(j));

				str(j)=convertCharsToStrings(chr(j));

				str(j)=strcat("x"+str(j));

		end

		X_base1 = convertStringsToChars(str(1));
		X_base2 = convertStringsToChars(str(2));
		

		TABLE2.Properties.RowNames={X_base1, X_base2,'cj-zj'};
		fprintf("\n\n");
		fprintf(" ----------------- Epomenos Simplex pinakas ----------------- \n\n");
		disp(TABLE2);
	else			

	% An ola ta stoixeia tis Cj-Zj grammhs den einai thetika, tote exei vrethei h veltisti lusi
		
		flag=false;
		
		fprintf('\nBrethike h veltisti lusi. Einai h z = %d, me times:\n\n', Solution(end));
		
		disp(Current_Solution);
	end
end
//Cálculo dos parâmetros da PC-SAFT usando recaracterização de Lage (2007), 5 pseudocomponentes, de Al Ajmi 2011: F1
clear
clc

exec("C:\Users\jvms9\OneDrive\Área de Trabalho\Comparação Al-Ajmi CBTERMO\Etapa 2\F1\Lage, 5 pseudos\Função PC-SAFT para obter Z.sce")

//Informações dos pseudocomponentes gerados na recaracterização de Lage (2007)
pmolar = [14.653409628;16.245857733;8.967827490;3.318413271;6.659491879] //%molar
xf = pmolar/100
xn = xf./sum(xf)
M = [107.796529830;185.413427650;314.719380540;448.625722000;841.210348370]//Massa molar dos pseudocomponentes

x7mais = sum(xf)
M7mais = 291
SG7mais = 0.8945/0.999
xn7mais = sum(xn)//fração molar 7+ normalizada

//Cálculo de SG de cada fração através da correlação de Soreide
function z = s(Cf)
    for i = 1:length(M)
        SG(i) = 0.2855 + Cf*(M(i)-66)^0.13
    end
    soma5 = 0
    for i = 1:length(M)
        soma5 = soma5 + xn(i)*M(i)/SG(i)
    end
    z = SG7mais - xn7mais*M7mais/soma5
endfunction
Cfchute = 0.5
Cf = fsolve(Cfchute,s)
SG = 0.2855 + Cf*(M-66)^0.13
//Adicionando SG e M do C6 de acordo com a tabela 4.6 de Riazi (2005) p.162
SG = [0.690;SG]
M = [84;M]

//Cálculo de propriedades como densidade, índice de refração de cada pseudocomponente
In = 0.34-exp(2.30884-2.96508*M^0.1) //Parâmetro do índice de refração p.161 - 162 Riazi (2005) 20°C
nr = sqrt((1+2*In)./(1-In)) //p. 66 Riazi (2005), índice de refração a 20°C
nr15 = nr - 0.0004*(273.15+15.5-293.15) //p. 67 Riazi (2005), índice de refração a 15,5°C
In15 = (nr15^2-1)./(nr15^2+2)
d15 = 0.999*SG //densidade a 15,5°C
d = In.*d15./In15 //densidade a 20°C p.68 Riazi (2005)

//Cálculo das frações PNA ref: p.126-127 Riazi (2005) (método ndm)
v = 2.51*(nr-1.475)-(d-0.851)
w = (d-0.851)-1.11*(nr-1.475)
for i = 1:length(SG)
    if v(i) > 0 then
        a(i) = 430
        b(i) = 0.055
    else
        a(i) = 670
        b(i) = 0.080
    end
end
S = zeros(length(v),1) //teor de enxofre
for i = 1:length(SG)
    if w(i) > 0 then
        pCr(i) = 820*w(i)-3*S(i)+10000/M(i) //percentual de carbono em aneis
    else
        pCr(i) = 1440*w(i)-3*S(i)+10600/M(i) //percentual de carbono em aneis
    end
end
pCa = a.*v + 3660./M //Percentual de carbono aromático
pCn = pCr - pCa //Percentual de carbono naftênico
pCp = 100 - pCr //Percentual de carbono parafínico
xCa = pCa/100 //Fração de carbono aromático
xCn = pCn/100 //Fração de carbono naftênico
xCp = pCp/100 //Fração de carbono parafínico

//Corrigindo frações negativas que ocorreram para SCN com M<200 (a fração negativa é descontada de todas as frações, após isso, é feita rernormalização)
PNA = [xCp,xCn,xCa]
[linh,col] = size(PNA)
for i = 1:linh
    for j = 1:col
        if PNA(i,j)<0 then
            const = PNA(i,j)
            PNA(i,1) = PNA(i,1)-const
            PNA(i,2) = PNA(i,2)-const
            PNA(i,3) = PNA(i,3)-const
            soma = sum(PNA(i,:))
            PNA(i,1) = PNA(i,1)/soma
            PNA(i,2) = PNA(i,2)/soma
            PNA(i,3) = PNA(i,3)/soma
        end
    end
end
xCp = PNA(:,1)
xCn = PNA(:,2)
xCa = PNA(:,3)

//Cálculo dos parâmetros m, epsilon_K e sigma da PC-SAFT para cada número de carbono (LIANG et al., 2014)
for i = 1:length(SG)
    m(i) = (0.02569*M(i)+0.8709)*xCp(i) + (0.02254*M(i)+0.6827)*xCn(i) + (0.02576*M(i)+0.2588)*xCp(i)
    me_K(i) = 6.8311*M(i)+124.42 //Eq.4 Liang et al. (2014)
    e_K(i) = me_K(i)/m(i)
end

//Cálculo do parâmetro sigma da PC-SAFT para cada pseudocomponente
T = 273.15+15.5 //K temperatura referente à densidade
P = 101325 //Pa
R = 8.3144621// J K^-1 mol^-1

//Cálculo de densidade a 15,5ºC e Vmolar
df = d15.*(1./(M*1e-6)) //convertendo para mol/m3
Vmolar = 1./df //m3/mol

K = 1.380649e-23 //m^2 kg s^-2 K^-1
Na = 6.02214076e23 //mol^-1

for i = 1:length(SG)
    x = 1
    sigmachute = 3.71 //angstrons, estimativa
    epsilon_K(i) = e_K(i) //K
    k = 0
    function z = g(sigma)
        ro(i) = (Vmolar(i)*1e10^3/Na)^-1//angstrons^-3
        d(i) = sigma*(1-0.12*exp(-3*epsilon_K(i)/T)) //angstrons
        eta(i) = %pi/6*ro(i)*x*m(i)*d(i)^3
        Z = SAFTZ(eta(i),T,x,m(i),sigma,epsilon_K(i),k)
        z = P-Z*R*T/Vmolar(i)
    endfunction
    sigma(i) = fsolve(sigmachute,g)
end


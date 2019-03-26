clc; close all; clear all;

load('C:\Users\Fernando\Dropbox\Mestrado\Simulacoes\simulações e testes\bases_rossler\rossler_2_8');
t = 0:deltat:(npontos-1)*deltat;

%GRÁFICOS

%Figura 2.8

%Figura A)

    %---------------------------------%
    %-- Atrator de Rössler - x vs y --%
    %---------------------------------%
    figure(1);
    plot(x(1:100000,1),x(1:100000,2))
    xlabel('x')
    ylabel('y')
    grid on;

%Figura B)

    %------------------------------------------------------%
    %-- log PDS - Espectro de frequência da trajetória x --%
    %------------------------------------------------------%
    figure(2);
    [pxx,w] = periodogram(x(:,1));%,window);
    plot(w,10*log10(pxx));
    xlabel('\omega')
    ylabel('log(PSD)')
    grid on;

%Figura C)

    %---------------------------------%
    %-- Evolução temporal das fases --%
    %---------------------------------%
    
    %Seção de Poincaré
        %-- Pontos em que o valor de y é aproximadamente 0 e x < 0--%
        figure(3);
        y = x(:,2);
        yacumulado = cumsum(y);
        [picos,locs] = findpeaks(yacumulado);
        subplot(211); plot(t,y,'k',t(locs), y(locs),'go');grid;%xlim([0 2000]);
        subplot(212); plot(t,yacumulado,'k',t(locs), yacumulado(locs),'go');grid;%xlim([0 2000]);
        
        %-- Tempos em que os pontos da figura anterior ocorrem --%
        figure(4);
        plot(x(:,1),x(:,2),'k',x(locs,1),x(locs,2),'go');
        tempospoinc = locs;
        
        figure(5);
        plot(locs)
        
        %-- Cálculo da fase pela seção de Poincaré --%
        
 % Será gerado um loc artificial
locs(331) = 100000;

   k = 1;
   j = 1;
   r = 1;
        
while j < length(locs)
    
   if t(k) < t(locs(j))
    
   fi_p(r) = 2*pi()*(t(k) - t(locs(j)))/((t(locs(j+1)))-t(locs(j)))+2*pi()*j;
   
   r = r+1;
   k = k+1;
   
   else
   
   j = j+1;
   
   end
   
end   

for r = 199281:200000
fi_p(r) = ((2*pi()*(t(k) - t(locs(330)))/((t(locs(330)))-t(locs(331)))+2*pi()*j));
end

phi_p = fi_p - 2*pi;
figure(6);
plot(phi_p)
xlabel('t')
ylabel('\fi_p')


    t = 1:1000;
    %Função cos(2*pi*t)
    %x = cos(2*pi*t)
    %Transformada de Hilbert 
       
    figure(7);
    xh = hilbert(x(:,1));
    fi_xh = unwrap(angle(xh)); 
    plot(t,unwrap(angle(xh)))
    
    %Arco tangente
    figure(8);
    fi = atan2(x(:,2),x(:,1));
    fi_arctan_acum = unwrap(fi);
    plot(t,fi,t,fi_arctan_acum)
    
    %Figura D)
    
        %-----------------------------------%
        %-- Gráfico das fases sobrepostas --%
        %-----------------------------------%
        
plot(t,fi_p,t,fi_arctan_acum,'-.',t,fi_xh,'--')

figure(9);
subplot(221)
plot(x(1:50000,1),x(1:50000,2));
xlabel('x');
ylabel('y');
grid on;

subplot(222)
plot(100*w,10*log10(pxx));
xlabel('\omega');
ylabel('log(PSD)');
xlim([0 5])
ylim([-30 70])
grid on;

subplot(223)
plot(t,fi_p,t,fi_arctan_acum,'-.',t,fi_xh,'--')
xlabel('t');
ylabel('\phi');
xlim([0 1000])
ylim([0 1000])
grid on;

subplot(224)
plot(t,fi_p,t,fi_arctan_acum,'-.',t,fi_xh,'--')
xlabel('t');
ylabel('\phi');
xlim([210 220])
ylim([210 240])
grid on;

% Tentativa de encontrar as frequências em que ocorre o sincronismo de fase

for freq = 1:length(locs)-1
   
    w(freq) = 2*pi/(t(locs(freq+1))-t(locs(freq)));
    
freq        
end

    figure(10)
    npontos = 100;
    t = 0:deltat:(npontos-1)*deltat;
    %Função cos(2*pi*t)
    x = cos(2*pi*t);
    subplot(211);
    plot(t,x)
    xlabel('t')
    ylabel('x(t)')
    %Transformada de Hilbert 
    subplot(212)
    xh = hilbert(x);
    fi_xh = unwrap(angle(xh)); 
    plot(t,unwrap(fi_xh)/pi)
    xlabel('t')
    ylabel('\frac{\phi}{\pi}')







    %Seção de Poincaré
        %-- Pontos em que o valor de y é aproximadamente 0 e x < 0--%
        figure(3);
        y = x(:,2);
        yacumulado = cumsum(y);
        [picos,locs] = findpeaks(yacumulado);
        plot(t,y,'k',t(locs), y(locs),'go');grid;%xlim([0 2000]);
        xlabel('t')
        ylabel('x(t)')
        %subplot(212); plot(t,yacumulado,'k',t(locs), yacumulado(locs),'go');grid;%xlim([0 2000]);


        figure(4);
        plot(x(1:10000,1),x(1:10000,2),'k',x(locs(1:20),1),x(locs(1:20),2),'go');
        xlabel('x(t)')
        ylabel('y(t)')
        tempospoinc = locs;
        
        figure(5);
        plot(locs(1:20))










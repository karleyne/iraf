###################################################################
# task original de Chico Jablonski - 2004
# 
# modificada por CVR (11/2004) para considerar barras de erros
# na entrada. Mas , o erro dos bins continua sendo calculado via dispersao. Especificamente,
# o erro de um bin eh o desvio padrao dos pontos no bin (NAO eh o desvio padrao da media).
#
# Karleyne: modificado em 08/2008 para converter mag em fluxo para a opção binada
#
###################################################################
procedure diagfase(arquivo)

char	arquivo	{prompt="Arquivo de dados"}
char	colunas	{prompt="Colunas com t,y e erro (se existir)"}
char	arqsai  {prompt="Nome do arquivo de saida (nulo p/ nao salvar)"}
real	periodo	{ 0.0605163312, prompt="Periodo"}
real	epoca	{ 2452969.322083, prompt="Epoca"}
bool	erros	{prompt="Considerar erros? y/n"}
real	Tshift	{0.,prompt="Constante a somar aos tempos"}
real	yshift	{0.,prompt="Constante a somar aos y"}
real	wx1	{0.,prompt="Limite inferior em fase"}
real	wx2	{2.,prompt="Limite superior em fase"}
real	xt	{0.5,prompt="Distancia para label tick em fase"}
real	xtp	{0.1,prompt="Distancia para tick pequeno em fase"}
real	wy1	{prompt="Limite inferior em y"}
real	wy2	{prompt="Limite superior em y"}
real	yt	{0.5, prompt="Distancia para label tick em y"}
real	ytp	{0.1, prompt="Distancia para tick pequeno em y"}
char	legx	{"Fase (ciclos)",prompt="Legenda x"}
char	legy	{prompt="Legenda y"}
char	titulo	{"",prompt="Titulo"}
char	cmd1	{"",prompt="Comando adicional do igi"}
char	cmd2	{"",prompt="Comando adicional do igi"}
char	simbolo	{"20 3", prompt="Simpolo para os pontos"}
char	lweight {"2",prompt="Espessura das linhas (1,2,3...)"}
char	arq2	{"",prompt="Arquivo a superpor"}
int	Nbins	{0,prompt="Se > zero, bina os dados"}
bool    conv    {yes, prompt="Coverte magnitude em Fluxo?"} 
real	f0	{3781., prompt="Flux of the zero magnitude"}
char	device	{"stdgraph",prompt="Device de saida"}


begin

char dados, dados_fase, dados_fase_sort,arqigi,lw
real x,y,err,tfase,fase,fase1,xmeio,ymeio,yabaixo,xaolado,dy,dfase,xlim,mmax,nmax
real sf[200],yf[200]
int i,ifase, nf[200],nb2
bool eerros

dados=mktemp('tmp')
dados_fase=mktemp('tmp')
dados_fase_sort=mktemp('tmp')
arqigi = mktemp('tmp') // ".igi"


eerros=erros
lw = lweight
if(device == 'stdgraph') lw="1"
fields(arquivo, fields=colunas, > dados)
#
# ######################
# zera somatorias
#
if (Nbins > 0){
   nb2=Nbins+Nbins
   for (i=1; i<=nb2; i=i+1){
      yf[i]=0.0
      sf[i]=0.0
      nf[i]=0
   }
}

list=dados

#
#se considera os erros entra neste if
#

if (erros) {
while (fscan(list,x,y,err) != EOF){
	if (Tshift != 0.0) x = x + Tshift
	if (yshift != 0.0) y = y + yshift
        if (conv){
          y = 10**(-y / 2.5) * f0
          err=0.921 * y * err 
        }
	tfase = (x-epoca)/periodo
	ifase = int(tfase)
	fase = tfase-ifase
	if (fase < 0.0) fase = fase + 1.0
	fase1 = fase+1.
	print(fase,' ',y,Ž Ž,err, >>dados_fase)
	if (fase1<wx2) print(fase1,' ',y,Ž Ž,err, >>dados_fase)
	if (Nbins > 0){
	   i = int(fase*Nbins+0.5) + 1
           yf[i]=yf[i]+y
	   sf[i]=sf[i]+y**2
	   nf[i]=nf[i]+1
   	   }
   	#print (x,y,err)
   	}
   	}

#
#se não considera erros vem por aqui
#
else {
while (fscan(list,x,y) != EOF){
	if (Tshift != 0.0) x = x + Tshift
	if (yshift != 0.0) y = y + yshift
        if (conv) {
           y = 10**(-y / 2.5)
        }
	tfase = (x-epoca)/periodo
	ifase = int(tfase)
	fase = tfase-ifase
	if (fase < 0.0) fase = fase + 1.0
	fase1 = fase+1.
	print(fase,' ',y, >>dados_fase)
	if (fase1<wx2) print(fase1,' ',y,' ', >>dados_fase)
	if (Nbins > 0){
	   i = int(fase*Nbins+0.5) + 1
           yf[i]=yf[i]+y
	   sf[i]=sf[i]+y**2
	   nf[i]=nf[i]+1
   	   }
   	#print (x,y,err)
        }
}

list=""

##########################
#se faz a binagem calcula as medias, sigmas
#

if (Nbins > 0){
   dfase = 1./Nbins
   yf[1]=yf[1]+yf[Nbins+1]
   sf[1]=sf[1]+sf[Nbins+1]
   nf[1]=nf[1]+nf[Nbins+1]
   for (i=1; i<=Nbins; i=i+1){
      if (nf[i]>0){
         yf[i]=yf[i]/nf[i]
         sf[i]=sqrt((sf[i]/nf[i]-yf[i]**2)/nf[i])
         if (nf[i]==1 ) sf[i]=yf[i]
         if (sf[i]==0.) sf[i]=yf[i]
         yf[Nbins+i]=yf[i]
         sf[Nbins+i]=sf[i]
         nf[Nbins+i]=nf[i]
      }     
   }

# ##############################
#
# normalizando o fluxo e o sigma fluxo
#
if(conv){
  mmax=0
  for (i=1; i<=Nbins;i=i+1)       mmax=max(mmax,yf[i])
  for (i=1; i<=nb2;i=i+1){ 
      yf[i]=yf[i]/mmax
      sf[i]=sf[i]/mmax
      print (i,yf[i],sf[i],nf[i])
      }
  }  
#
###################################
# 
#  Criando arquivo de saida 
#
   printf("%10.5f %10.5f %10.5f\n",0.0,yf[1],sf[1], >dados_fase_sort)
   #print(0.,' ',yf[1],' ',sf[1], >dados_fase_sort)
   for(i=2; i<=nb2; i=i+1){
      fase = (i-1)*dfase
      printf("%10.5f %10.5f %10.5f\n",fase,yf[i],sf[i], >>dados_fase_sort)
      #print(fase,' ',yf[i],' ',sf[i], >>dados_fase_sort)  
   }
}
else{  
   sort(dados_fase,col=1,num+, >dados_fase_sort)
}

if (arqsai != "") copy(dados_fase_sort,arqsai,ver-)

if (Nbins > 0) eerros = yes

print('erase', >arqigi)
#print('vpage 0.15 0.85 0.15 0.85', >>arqigi)
print('lwei ',lw, >>arqigi)
print('lim ',wx1,' ',wx2,' ',wy1,' ',wy2, >>arqigi)
print('tick ',xtp,' ',xt,' ',ytp,' ',yt, >>arqigi)
print('pt ',simbolo,' ; exp 0.25', >>arqigi)
print('data ',dados_fase_sort, >>arqigi)
print('xcol 1; ycol 2', >>arqigi)
if (eerros){
      print('ecol 3', >>arqigi)
      }
#if (Nbins > 0){
#      print('ecol 3', >>arqigi)
#      }
print('point ', >>arqigi)
if (eerros){
   print('errorb -2 ; errorb 2', >>arqigi)
}
#if (Nbins > 0){
#   print('errorb -2 ; errorb 2', >>arqigi)
#}
if (access(arq2)){
   print('data ',arq2, >>arqigi)
   print('xcol 1 ; ycol 2', >>arqigi)
   print('connect ', >>arqigi)
}
print('exp 1.0 ; box', >>arqigi)
xmeio = (wx2+wx1)/2.
ymeio = (wy2+wy1)/2.
dy = (wy2-wy1)/10.
yabaixo = wy1 - dy*0.8
xaolado = wx1 - (wx2-wx1)/11.
print('reloc ',xmeio,' ',yabaixo,'put 5 \\r',legx,'\e', >>arqigi)
print('angle 90.', >>arqigi)
print('reloc ',xaolado,' ',(ymeio+dy*0.0),'put 5 \\r',legy,'\e', >>arqigi)
print('angle 0.', >>arqigi)
if (cmd1 != "") print(cmd1, >>arqigi)
if (cmd2 != "") print(cmd2, >>arqigi)
print('exp 1.2 ; title ',titulo, >>arqigi)
print('end', >>arqigi)

type(arqigi) | igi(device=device)

gflush

delete(dados,ver-)
delete(dados_fase,ver-)
delete(dados_fase_sort,ver-)
delete(arqigi,ver-)

end

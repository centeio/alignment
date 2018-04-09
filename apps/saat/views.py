"""Views for the base app"""

from django.shortcuts import render
from .models import SeqForm
from django.shortcuts import redirect
from .algorithms.AffineGlobal import AffineGlobal
from .algorithms.AffineLocal import AffineLocal
from .algorithms.Ebi import Ebi
from .algorithms.util.alignment import Parameters


def local(request):
    if request.method == "POST":
        form = SeqForm(request.POST)
        
        if form.is_valid():
            savedForm = form.save(commit=False)

            algorithm = AffineLocal(multiples=True)

            ebi = Ebi(type="LOCAL")

            gapopen = savedForm.gapopen
            gapext = savedForm.gapext

            #Definicoes dos parametros
            par = Parameters(gapopen=-gapopen,matrix=savedForm.matrix,stype=savedForm.seqtype)

            #Sequencias
            seq1=savedForm.seq1.upper()
            seq2=savedForm.seq2.upper()
            
            result_from_algo = algorithm.run(seq1, seq2, par)

            #EBI
            result_from_ebi = ebi.run(seq1, seq2, savedForm.seqtype, savedForm.matrix, gapopen, gapext)
            
            return render(request, 'base/result.html', {'algo_results': result_from_algo, 'ebi': result_from_ebi})
    else:
        form = SeqForm()
        
    return render(request, 'base/home.html', {'form': form, 'type':'Local Alignment with Affine Gap Penalty'})

def home(request):
    if request.method == "POST":
        form = SeqForm(request.POST)
        if form.is_valid():
            savedForm = form.save(commit=False)

            algorithm = AffineGlobal(multiples=True)
            ebi = Ebi(type="GLOBAL")

            gapopen = savedForm.gapopen
            gapext = savedForm.gapext


            #Definicoes dos parametros
            par = Parameters(gapopen=-gapopen,gapext=-gapext,matrix=savedForm.matrix,stype=savedForm.seqtype)

            #Sequencias
            seq1=savedForm.seq1.upper()
            seq2=savedForm.seq2.upper()
            
            result_from_algo = algorithm.run(seq1, seq2, par)

            #EBI
            result_from_ebi = ebi.run(seq1, seq2, savedForm.seqtype, savedForm.matrix, gapopen, gapext)
            
            return render(request, 'base/result.html', {'algo_results': result_from_algo, 'ebi': result_from_ebi})
    else:
        form = SeqForm()
        
    return render(request, 'base/home.html', {'form': form, 'type':'Global Alignment with Affine Gap Penalty'})

def result(request):
    return render(request, 'base/result.html', {'result': result})

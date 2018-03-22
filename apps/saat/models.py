"""Base models"""
from django.db import models
from django.utils import timezone
from django import forms



MATRICES = (
    (1, ("BLOSUM62")),
    (2, ("PAM250")),
    (3, ("DNAFULL"))
)

SEQ_TYPE = (
    (1, ("PROTEIN")),
    (2, ("DNA"))
)

GAP_VALUES = (
    (1, (1)),
    (5, (5)),
    (10, (10)),
    (15, (15)),
    (20, (15)),
    (50, (50)),
    (100, (100))  
)


GAP_EXT_VALUES = (
    (0.0005, (0.0005)),
    (0.001, (0.001)),
    (0.05, (0.05)),
    (0.1, (0.1)),
    (0.2, (0.2)),
    (0.4, (0.4)),
    (0.5, (0.5)),
    (0.6, (0.6)),
    (0.8, (0.8)),
    (1.0, (1.0)),
    (5.0, (5.0)),
    (10.0, (10.0))  
)


class Alignment(models.Model):
    seq1 = models.TextField()
    seq2 = models.TextField()
    matrix = models.IntegerField(choices=MATRICES, default=1)
    seqtype = models.IntegerField(choices=SEQ_TYPE, default=1)
    gapopen = models.IntegerField(choices=GAP_VALUES, default=10)
    gapext = models.DecimalField(max_digits=3, decimal_places=1, choices=GAP_EXT_VALUES, default=0.5)

    
    def __str__(self):
        return "seq1: {}\n seq2:{}\n matrix:{}\n seqtype:{}\n gapopen:{}\n gapext:{}".format(
            self.seq1, self.seq2, self.matrix, self.seqtype, self.gapopen, self.gapext
            )

class SeqForm(forms.ModelForm):

    class Meta:
        model = Alignment
        fields = ('seq1', 'seq2', 'seqtype', 'matrix', 'gapopen', 'gapext')
As usual there is no one unique solution ! 

Here is an example:

```{.python}
def aln_Convert(input,  output, input_format="clustal", output_format="stockholm"):
   alignments = AlignIO.parse(input, input_format)
   count = AlignIO.write(alignments, output, output_format)
   print("Converted %i alignments" % count)


aln_Convert("data/muscle-patato_pep.clw","muscle-patato_pep.stk")

os.system("head -4 muscle-patato_pep.stk")

```


> \# STOCKHOLM 1.0   
> \#=GF SQ 10   
> PGSC0003DMP400022612 -------------------------------------------------------------------------------------------------------------------------------MLPFESIEEASMSLGRNLTFGETLWF----------------------------NYTADKSDFYLYCHNTIFLI--------------------------------------------------------------------------------------------------------------------------------------IF------------------YSLV--PLPMVMIEL-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------LRSKKFEI-------------------------------------------------------------------------YKIQPKTTLTLSEMME-----------------------CYRKVLMTFVFAVG-------P------LQLFSYPIIEWVGIRTSLPLPSASEVFWQLVVYFV---------------------------------------------------------------------------------------VEDYANYWLHRMLHHYK--------------------------------------------------------------------WGY--DKIHSVHHEYVAPISFAAPYAHWAE-----IIILGLASFLGPLLVPCHMFTFLLWFVLR--------QIEAIETHSGYEFPWSPS-------------------KYIPFY-----------------------------------------------------------GGAIYHDY---HHFVGESCHSNFASVFTYCDYIYGTDKGFRYQKEVFEKRREKLRMEE---------------KVNGSAPLFKSE-----------   
> \#=GS PGSC0003DMP400022612 AC PGSC0003DMP400022612   



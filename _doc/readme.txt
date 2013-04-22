Ce template a �t� r�alis� � partir de la classe THLORIA.cls.
N'h�sitez donc pas � regarder la doc TL-user.pdf, c'est tr�s bien expliqu�!


Voil� quelques compl�ments.
Si vous avez la moindre question, n'h�sitez pas � me la poser : cedric.roux@gmail.com
Bonne r�daction!
(et si vous d�cidez d'utiliser ce template, n'h�sitez pas � me le dire, �a me fera plaisir!)


IMPORTANT
*********
- Les titres des chapitres sont � modifier dans le fichier mep/macros.tex

- pour la liste des notations, il faut utiliser makeindex avec les arguments suivants :
	"%bm".nlo -s nomencl.ist -o "%bm".nls
  Dans Texniccenter sous Windows, il suffit d'aller dans Build/Define Output Profiles, et copier/coller la ligne pr�c�dente dans la case "Command line arguments to passe to MakeIndex" (et du coup compiler N fois avec N tr�����s grand devant 1)


MOINS IMPORTANT
***************
- pour le fichier biblio, je vous conseille d'utiliser JabRef sous windows, un bon �diteur de fichier bib.

- FaitMinitocs : pour avoir des minitables des mati�res au d�but de chaque chapitre (1) ou pas (0)

- Figures : \bfigsh et \efigsh pour une figure ombr�e, \bfig et \efig pour une figure encadr�e (si \CadresOvales=1 dans le fichier these.tex), \bfigp et \efigp pour une figure pleine page

- FaitNomenclature : pour avoir la liste des notations utilis�es ou non. La syntaxe est \nome{$notation$}{Description}

- IncludePart... : pour gagner du temps de compilation!

- De nombreuses bo�tes sont d�finies (en plus de \AN{}, \Remarque{} et \Resultat{} qui sont illustr�es ici), voir le fichier macros pour plus de d�tails

- \WriteThisInToc\listoffigures : le \WriteThisInToc permet d'�crire "liste des figures" dans la table des mati�res

- Les commandes \AjouteLigne et \RetireLigne sont tr�s utiles!!! (cf fichier macro pour plus de d�tails)
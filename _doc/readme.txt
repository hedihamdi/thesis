Ce template a été réalisé à partir de la classe THLORIA.cls.
N'hésitez donc pas à regarder la doc TL-user.pdf, c'est très bien expliqué!


Voilà quelques compléments.
Si vous avez la moindre question, n'hésitez pas à me la poser : cedric.roux@gmail.com
Bonne rédaction!
(et si vous décidez d'utiliser ce template, n'hésitez pas à me le dire, ça me fera plaisir!)


IMPORTANT
*********
- Les titres des chapitres sont à modifier dans le fichier mep/macros.tex

- pour la liste des notations, il faut utiliser makeindex avec les arguments suivants :
	"%bm".nlo -s nomencl.ist -o "%bm".nls
  Dans Texniccenter sous Windows, il suffit d'aller dans Build/Define Output Profiles, et copier/coller la ligne précédente dans la case "Command line arguments to passe to MakeIndex" (et du coup compiler N fois avec N trèèèèès grand devant 1)


MOINS IMPORTANT
***************
- pour le fichier biblio, je vous conseille d'utiliser JabRef sous windows, un bon éditeur de fichier bib.

- FaitMinitocs : pour avoir des minitables des matières au début de chaque chapitre (1) ou pas (0)

- Figures : \bfigsh et \efigsh pour une figure ombrée, \bfig et \efig pour une figure encadrée (si \CadresOvales=1 dans le fichier these.tex), \bfigp et \efigp pour une figure pleine page

- FaitNomenclature : pour avoir la liste des notations utilisées ou non. La syntaxe est \nome{$notation$}{Description}

- IncludePart... : pour gagner du temps de compilation!

- De nombreuses boîtes sont définies (en plus de \AN{}, \Remarque{} et \Resultat{} qui sont illustrées ici), voir le fichier macros pour plus de détails

- \WriteThisInToc\listoffigures : le \WriteThisInToc permet d'écrire "liste des figures" dans la table des matières

- Les commandes \AjouteLigne et \RetireLigne sont très utiles!!! (cf fichier macro pour plus de détails)
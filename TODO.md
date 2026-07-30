# TODO

## Statistiques prévues

### sum-of-pairs (`-s sumofpairs`)

Statistique par paire de colonnes, basée sur `Stat2D` (déjà présente dans
`statistic.h`, actuellement sans implémentation concrète). Retirée de
l'aide (`--help`) tant qu'elle n'est pas implémentée, pour ne pas laisser
un utilisateur taper `-s sumofpairs` et tomber sur une erreur "not found".

À faire pour l'intégrer :
- Créer `sumofpairs.h`/`sumofpairs.cpp`, sur le modèle des 6 statistiques
  existantes (`gap`, `kabat`, `wentropy`, `trident`, `jensen`, `mvector`).
- Hériter de `Stat2D` plutôt que `Stat1D` (résultat par paire de colonnes,
  pas par colonne seule).
- Enregistrer la nouvelle classe dans `AddAllStatistics()` (`statistic.cpp`).
- Ajouter `sumofpairs` à la liste des statistiques dans `print_usage()`
  (`options.h`).
- Ajouter `tests/test_sumofpairs.cpp`, sur le modèle des fichiers de test
  existants (fixture synthétique, valeurs nominales vérifiées à la main
  ou par une implémentation Python indépendante, mode `--global`).

private static void ApplySynergy(PlayerCharacter ch, RaceClassSynergy synergy)
{
    ch.Strength     += synergy.BonusStrength;
    ch.Agility      += synergy.BonusAgility;
    ch.Intelligence += synergy.BonusIntelligence;
}

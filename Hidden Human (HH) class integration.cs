public static bool IsHiddenHumanUnlocked(PlayerProfile profile, LineageData lineage)
{
    // Example conditions:
    return lineage.IsJarvondis || profile.HasCompletedQuest("SEE_THE_UNSEEN");
}

[Serializable]
public class PlayerCharacterSave
{
    public string Name;
    public Race Race;
    public ClassType Class;
    public int Strength;
    public int Agility;
    public int Intelligence;

    public AppearanceData Appearance;
    public SkillNode[] SkillTree;
    public BackstoryData Backstory;
    public LineageData Lineage;
}

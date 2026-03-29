public enum Race { Human, Elf, Orc, HiddenHuman /* HH */, ... }
public enum ClassType { Warrior, Mage, ShadowPanther, ... }

public class CharacterTemplate
{
    public string Name;
    public Race Race;
    public ClassType Class;
    public int Strength;
    public int Agility;
    public int Intelligence;

    public SkillNode[] SkillTreeTemplate;
    public AppearanceTemplate DefaultAppearance;
    public LineageData LineageTemplate;
}

public class PlayerCharacter : CharacterTemplate
{
    public bool IsCustom;
    public AppearanceData Appearance;
    public SkillNode[] SkillTree;   // Instanced from template
    public BackstoryData Backstory;
}

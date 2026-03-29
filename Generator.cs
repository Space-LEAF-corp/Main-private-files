public static class BackstoryGenerator
{
    public static BackstoryData Generate(Race race, ClassType cls, LineageData lineage)
    {
        var data = new BackstoryData();

        // Very simple example:
        if (lineage.IsJarvondis)
        {
            data.Title = "Heir of Jarvondis";
            data.Summary = "Born under a hidden constellation, bound to the Panther’s shadow.";
        }
        else
        {
            data.Title = $"{race} {cls}";
            data.Summary = "A wanderer whose story is still being written.";
        }

        // Add key events based on race/class/lineage
        return data;
    }
}

---
title: "Main"
permalink: /posts/
author_profile: true
---

<div style="text-align: justify;">
This tool is a further refinement and development of lineage tracing. Lineage tracing refers to reconstructing cellular dynamics during the clonal process by inserting DNA base sequences. Step 1-3: Capture cell samples and insert fragments using DOX technology. Step 4-6: When the cells begin to clone, these fragments are carried along. As the cells differentiate, these fragments are also retained but may acquire some random mutations. This technology enables us to trace and observe cells. However, it is extremely expensive and has a very low success rate in DNA insert and detect. Therefore, using RNA sequencing for cell tracking has become a popular topic.
</div>

![图片](/assets/images/Main1.png)
<div style="text-align: justify;">
Due to the randomness of DNA fragment capture, this technology is difficult to apply to cells that have temporarily or permanently lost their differentiation capacity. Furthermore, existing algorithms cannot be used to distinguish the fate of mature cells, as their transcription profiles are too similar, and the accuracy is low even in datasets with clear differentiation boundaries. However, research related to mature cells is emerging and requires tracing technology. For example, ferroptosis may offer new methods for cancer treatment, but its mechanism remains unclear.
</div>

![图片](/assets/images/Main2.png)
You'll find this post in your `_posts` directory. Go ahead and edit it and re-build the site to see your changes. You can rebuild the site in many different ways, but the most common way is to run `jekyll serve`, which launches a web server and auto-regenerates your site when a file is updated.

To add new posts, simply add a file in the `_posts` directory that follows the convention `YYYY-MM-DD-name-of-post.ext` and includes the necessary front matter. Take a look at the source for this post to get an idea about how it works.

Jekyll also offers powerful support for code snippets:

```ruby
def print_hi(name)
  puts "Hi, #{name}"
end
print_hi('Tom')
#=> prints 'Hi, Tom' to STDOUT.
```

Check out the [Jekyll docs][jekyll-docs] for more info on how to get the most out of Jekyll. File all bugs/feature requests at [Jekyll’s GitHub repo][jekyll-gh]. If you have questions, you can ask them on [Jekyll Talk][jekyll-talk].

[jekyll-docs]: https://jekyllrb.com/docs/home
[jekyll-gh]:   https://github.com/jekyll/jekyll
[jekyll-talk]: https://talk.jekyllrb.com/
